import logging
import json
import os
import zipfile
from io import StringIO
import uuid
import shlex
import shutil
from datetime import datetime
from wsgiref.util import FileWrapper
from dateutil.parser import parse
import pytz
from pathlib import Path

import pandas as pd

from django.db import connections
from django.http import HttpResponse, FileResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.core.mail import send_mail
from django.conf import settings
from django.http import JsonResponse
from django.views import View

from django.urls import reverse

from rest_framework import status, viewsets, permissions
from rest_framework.exceptions import ParseError
from rest_framework.parsers import BaseParser
from rest_framework.response import Response
from rest_framework.views import APIView

from celery.result import AsyncResult

from api.security import ISpyBSafeQuerySet

from api.utils import get_params, get_highlighted_diffs, pretty_request
from viewer.utils import create_squonk_job_request_url
from viewer.utils import handle_uploaded_file

from viewer import filters
from viewer import models
from viewer import serializers
from viewer.squonk2_agent import Squonk2AgentRv, Squonk2Agent, get_squonk2_agent
from viewer.squonk2_agent import AccessParams, CommonParams, SendParams, RunJobParams



from .forms import CSetForm
from .tasks import (
    erase_compound_set_job_material,
    process_compound_set,
    process_design_sets,
    process_job_file_transfer,
    process_compound_set_job_file,
    validate_compound_set,
    task_load_target,
)
from .discourse import create_discourse_post, list_discourse_posts_for_topic, check_discourse_user
from .download_structures import (
    check_download_links,
    recreate_static_file,
    maintain_download_links
)

from .squonk_job_file_transfer import (
    check_file_transfer
)

from .squonk_job_request import (
    check_squonk_active,
    get_squonk_job_config,
    create_squonk_job,
)

logger = logging.getLogger(__name__)

# Fields injected in a session object to pass
# messages between views. This is used by UploadCSet
# to pass errors and other messages back to the user
# via the upload-cset.html template.
_SESSION_ERROR = 'session_error'
_SESSION_MESSAGE = 'session_message'

_SQ2A: Squonk2Agent = get_squonk2_agent()


class CompoundIdentifierTypeView(viewsets.ModelViewSet):
    queryset = models.CompoundIdentifierType.objects.all()
    serializer_class = serializers.CompoundIdentifierTypeSerializer
    permission_classes = [permissions.IsAuthenticated]


class CompoundIdentifierView(viewsets.ModelViewSet):
    queryset = models.CompoundIdentifier.objects.all()
    serializer_class = serializers.CompoundIdentifierSerializer
    permission_classes = [permissions.IsAuthenticated]
    filterset_fields = ["type", "compound"]


class VectorsView(ISpyBSafeQuerySet):
    """Vectors (api/vector)
    """
    queryset = models.SiteObservation.objects.all()
    serializer_class = serializers.VectorsSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class GraphView(ISpyBSafeQuerySet):
    """Graph (api/graph)
    """
    queryset = models.SiteObservation.objects.all()
    serializer_class = serializers.GraphSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class MolImageView(ISpyBSafeQuerySet):
    """Molecule images (api/molimg)
    """
    queryset = models.SiteObservation.objects.filter()
    serializer_class = serializers.MolImageSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filterset_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class CompoundImageView(ISpyBSafeQuerySet):
    """Compound images (api/cmpdimg)
    """
    queryset = models.Compound.objects.filter()
    serializer_class = serializers.CmpdImageSerializer
    filter_permissions = "project_id"
    filterset_fields = ("smiles",)


class ProteinMapInfoView(ISpyBSafeQuerySet):
    """Protein map info (file) (api/protmap)
    """
    queryset = models.SiteObservation.objects.filter()
    serializer_class = serializers.ProtMapInfoSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


class ProteinPDBInfoView(ISpyBSafeQuerySet):
    """Protein apo pdb info (file) (api/protpdb)
    """
    queryset = models.SiteObservation.objects.filter()
    serializer_class = serializers.ProtPDBInfoSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


class ProteinPDBBoundInfoView(ISpyBSafeQuerySet):
    """Protein bound pdb info (file) (api/protpdbbound)
    """
    queryset = models.SiteObservation.objects.filter()
    serializer_class = serializers.ProtPDBBoundInfoSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


class ProjectView(ISpyBSafeQuerySet):
    """Projects (api/project)
    """
    queryset = models.Project.objects.filter()
    serializer_class = serializers.ProjectSerializer
    # Special case - Project filter permissions is blank.
    filter_permissions = ""


class TargetView(ISpyBSafeQuerySet):
    """Targets (api/targets)
    """
    queryset = models.Target.objects.filter()
    serializer_class = serializers.TargetSerializer
    filter_permissions = "project_id"
    filterset_fields = ("title",)


class CompoundView(ISpyBSafeQuerySet):
    """Compounds (api/compounds)
    """
    queryset = models.Compound.objects.filter()
    serializer_class = serializers.CompoundSerializer
    filter_permissions = "project_id"
    filterset_fields = ("smiles", "current_identifier", "inchi", "long_inchi")


class SiteObservationView(ISpyBSafeQuerySet):
    queryset = models.SiteObservation.objects.filter()
    serializer_class = serializers.SiteObservationReadSerializer
    filter_permissions = "target_id__project_id"
    filterset_fields = ("code", "target_id", "target_id__title", "prot_type")


def react(request):
    """We "START HERE". This is the first API call that the front-end calls.
    """

    discourse_api_key = settings.DISCOURSE_API_KEY

    context = {}

    # Is the Squonk2 Agent configured?
    logger.info("Checking whether Squonk2 is configured...")
    sq2_rv = _SQ2A.configured()
    if sq2_rv.success:
        logger.info("Squonk2 is configured")
        context['squonk_available'] = 'true'
    else:
        logger.info("Squonk2 is NOT configured")
        context['squonk_available'] = 'false'

    if discourse_api_key:
        context['discourse_available'] = 'true'
    else:
        context['discourse_available'] = 'false'

    user = request.user
    if user.is_authenticated:
        context['discourse_host'] = ''
        context['user_present_on_discourse'] = 'false'
        # If user is authenticated and a discourse api key is available, then check discourse to
        # see if user is set up and set up flag in context.
        if discourse_api_key:
            context['discourse_host'] = settings.DISCOURSE_HOST
            error, error_message, user_id = check_discourse_user(user)
            if user_id:
                context['user_present_on_discourse'] = 'true'

        # If user is authenticated Squonk can be called then return the Squonk host
        # so the Frontend can navigate to it
        context['squonk_ui_url'] = ''
        if sq2_rv.success and check_squonk_active(request):
            context['squonk_ui_url'] = _SQ2A.get_ui_url()

    return render(request, "viewer/react_temp.html", context)


def save_pdb_zip(pdb_file):
    zf = zipfile.ZipFile(pdb_file)
    zip_lst = zf.namelist()
    zfile = {}
    zfile_hashvals = {}
    print(zip_lst)
    for filename in zip_lst:
        # only handle pdb files
        if filename.split('.')[-1] == 'pdb':
            f = filename.split('/')[0]
            save_path = os.path.join(settings.MEDIA_ROOT, 'tmp', f)
            if default_storage.exists(f):
                rand_str = uuid.uuid4().hex
                pdb_path = default_storage.save(save_path.replace('.pdb', f'-{rand_str}.pdb'), ContentFile(zf.read(filename)))
            # Test if Protein object already exists
            # code = filename.split('/')[-1].replace('.pdb', '')
            # test_pdb_code = filename.split('/')[-1].replace('.pdb', '')
            # test_prot_objs = Protein.objects.filter(code=test_pdb_code)
            #
            # if len(test_prot_objs) > 0:
            #     # make a unique pdb code as not to overwrite existing object
            #     rand_str = uuid.uuid4().hex
            #     test_pdb_code = f'{code}#{rand_str}'
            #     zfile_hashvals[code] = rand_str
            #
            # fn = test_pdb_code + '.pdb'
            #
            # pdb_path = default_storage.save('tmp/' + fn,
            #                                 ContentFile(zf.read(filename)))
            else:
                pdb_path = default_storage.save(save_path, ContentFile(zf.read(filename)))
            test_pdb_code = pdb_path.split('/')[-1].replace('.pdb', '')
            zfile[test_pdb_code] = pdb_path

    # Close the zip file
    if zf:
        zf.close()

    return zfile, zfile_hashvals


def save_tmp_file(myfile):
    """ Save file in temporary location for validation/upload processing
    """

    name = myfile.name
    path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
    tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

    return tmp_file


class UploadCSet(APIView):
    """Render and control viewer/upload-cset.html  - a page allowing upload of computed sets. Validation and
    upload tasks are defined in `viewer.compound_set_upload`, `viewer.sdf_check` and `viewer.tasks` and the task
    response handling is done by `viewer.views.ValidateTaskView` and `viewer.views.UploadTaskView`
    """

    def get(self, request):

        tag = '+ UploadCSet GET'
        logger.info('%s', pretty_request(request, tag=tag))
        logger.info('User=%s', str(request.user))
#        logger.info('Auth=%s', str(request.auth))

        # Any messages passed to us via the session?
        # Maybe from a redirect?
        # It so take them and remove them.
        session_error = None
        if _SESSION_ERROR in request.session:
            session_error = request.session[_SESSION_ERROR]
            del request.session[_SESSION_ERROR]
        session_message = None
        if _SESSION_MESSAGE in request.session:
            session_message = request.session[_SESSION_MESSAGE]
            del request.session[_SESSION_MESSAGE]

        # Only authenticated users can upload files
        # - this can be switched off in settings.py.
        user = self.request.user
        if not user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            context = {}
            context['error_message'] \
                = 'Only authenticated users can upload files' \
                  ' - please navigate to landing page and Login'
            return render(request, 'viewer/upload-cset.html', context)

        form = CSetForm()
        existing_sets = models.ComputedSet.objects.all()
        context = {'form': form,
                   'sets': existing_sets,
                   _SESSION_ERROR: session_error,
                   _SESSION_MESSAGE: session_message}
        return render(request, 'viewer/upload-cset.html', context)

    def post(self, request):

        tag = '+ UploadCSet POST'
        logger.info('%s', pretty_request(request, tag=tag))
        logger.info('User=%s', str(request.user))
#        logger.info('Auth=%s', str(request.auth))

        # Only authenticated users can upload files
        # - this can be switched off in settings.py.
        user = self.request.user
        if not user.is_authenticated and settings.AUTHENTICATE_UPLOAD:
            context = {}
            context['error_message'] \
                = 'Only authenticated users can upload files' \
                  ' - please navigate to landing page and Login'
            return render(request, 'viewer/upload-cset.html', context)

        form = CSetForm(request.POST, request.FILES)
        if form.is_valid():

            # Get all the variables needed from the form.
            # The fields we use will be based on the 'submit_choice',
            # expected to be one of V (validate), U (upload) or D delete
            choice = request.POST['submit_choice']

            # Generate run-time error if the required form fields
            # are not set based on the choice made...

            # The 'sdf_file' anf 'target_name' are only required for upload/update
            sdf_file = request.FILES.get('sdf_file')
            target = request.POST.get('target_name')
            update_set = request.POST.get('update_set')

            logger.info('+ UploadCSet POST choice="%s" target="%s" update_set="%s"', choice, target, update_set)

            # If a set is named the ComputedSet cannot be 'Anonymous'
            # and the user has to be the owner.
            selected_set = None
            if update_set and update_set != 'None':
                computed_set_query = models.ComputedSet.objects.filter(unique_name=update_set)
                if computed_set_query:
                    selected_set = computed_set_query[0]
                else:
                    request.session[_SESSION_ERROR] = \
                        'The set could not be found'
                    logger.warning('- UploadCSet POST error_msg="%s"', request.session[_SESSION_ERROR])
                    return redirect('upload_cset')

            # If validating or uploading we need a Target and SDF file.
            # If updating or deleting we need an update set (that's not 'None')
            if choice in ['V', 'U']:
                if sdf_file is None or target is None:
                    request.session[_SESSION_ERROR] = \
                        'To Validate or Upload' \
                        ' you must provide a Target and SDF file'
                    logger.warning('- UploadCSet POST error_msg="%s"', request.session[_SESSION_ERROR])
                    return redirect('upload_cset')
            elif choice in ['D']:
                if update_set == 'None':
                    request.session[_SESSION_ERROR] = \
                        'To Delete you must select an existing set'
                    logger.warning('- UploadCSet POST error_msg="%s"', request.session[_SESSION_ERROR])
                    return redirect('upload_cset')

            # If uploading (updating) or deleting
            # the set owner cannot be anonymous
            # and the user needs to be the owner
            if choice in ['U', 'D'] and selected_set:
                if selected_set.owner_user.id == settings.ANONYMOUS_USER:
                    request.session[_SESSION_ERROR] = \
                        'You cannot Update or Delete Anonymous sets'
                elif selected_set.owner_user != user:
                    request.session[_SESSION_ERROR] = \
                        'You can only Update or Delete sets you own'
                # Something wrong?
                # If so redirect...
                if _SESSION_ERROR in request.session:
                    logger.warning('- UploadCSet POST error_msg="%s"', request.session[_SESSION_ERROR])
                    return redirect('upload_cset')

            # Save uploaded sdf and zip to tmp storage
            tmp_pdb_file = None
            tmp_sdf_file = None
            if 'pdb_zip' in list(request.FILES.keys()):
                pdb_file = request.FILES['pdb_zip']
                tmp_pdb_file = save_tmp_file(pdb_file)
            if sdf_file:
                tmp_sdf_file = save_tmp_file(sdf_file)

            if choice == 'V':
                # Validate
                # Start celery task
                task_params = {'user_id': user.id,
                               'sdf_file': tmp_sdf_file,
                               'target': target}
                if tmp_pdb_file:
                    task_params['zfile'] = tmp_pdb_file
                if update_set:
                    task_params['update'] = update_set
                task_validate = validate_compound_set.delay(task_params)

                logger.info('+ UploadCSet POST "Validate" task underway')

                # Update client side with task id and status
                context = {'validate_task_id': task_validate.id,
                           'validate_task_status': task_validate.status}
                return render(request, 'viewer/upload-cset.html', context)

            elif choice == 'U':
                # Upload
                # Start chained celery tasks. NB first function passes tuple
                # to second function - see tasks.py
                task_params = {'user_id': user.id,
                               'sdf_file': tmp_sdf_file,
                               'target': target}
                if tmp_pdb_file:
                    task_params['zfile'] = tmp_pdb_file
                if update_set:
                    task_params['update'] = update_set
                task_upload = (
                        validate_compound_set.s(task_params) |
                        process_compound_set.s()).apply_async()

                logger.info('+ UploadCSet POST "Upload" task underway')

                # Update client side with task id and status
                context = {'upload_task_id': task_upload.id,
                           'upload_task_status': task_upload.status}
                return render(request, 'viewer/upload-cset.html', context)

            elif choice == 'D':
                # Delete
                selected_set.delete()

                request.session[_SESSION_MESSAGE] = \
                    f'Compound set "{selected_set.unique_name}" deleted'

                logger.info('+ UploadCSet POST "Delete" done')

                return redirect('upload_cset')

            else:
                logger.warning('+ UploadCSet POST unsupported submit_choice value (%s)', choice)

        else:
            logger.warning('- UploadCSet POST form.is_valid() returned False')

        logger.info('- UploadCSet POST (leaving)')

        context = {'form': form}
        return render(request, 'viewer/upload-cset.html', context)



def email_task_completion(contact_email, message_type, target_name, target_path=None, task_id=None):
    """Notify user of upload completion
    """

    logger.info('+ email_notify_task_completion: ' + message_type + ' ' + target_name)
    email_from = settings.EMAIL_HOST_USER

    if contact_email == '' or not email_from:
        # Only send email if configured.
        return

    if message_type == 'upload-success':
        subject = 'Fragalysis: Target: '+target_name+' Uploaded'
        message = 'The upload of your target data is complete. Your target is available at: ' \
                  + str(target_path)
    elif message_type == 'validate-success':
        subject = 'Fragalysis: Target: '+target_name+' Validation'
        message = 'Your data was validated. It can now be uploaded using the upload option.'
    else:
        # Validation failure
        subject = 'Fragalysis: Target: ' + target_name + ' Validation/Upload Failed'
        message = 'The validation/upload of your target data did not complete successfully. ' \
                  'Please navigate the following link to check the errors: validate_task/' + str(task_id)

    recipient_list = [contact_email, ]
    logger.info('+ email_notify_task_completion email_from: %s', email_from )
    logger.info('+ email_notify_task_completion subject: %s',  subject )
    logger.info('+ email_notify_task_completion message: %s',  message )
    logger.info('+ email_notify_task_completion contact_email: %s', contact_email )

    # Send email - this should not prevent returning to the screen in the case of error.
    send_mail(subject, message, email_from, recipient_list, fail_silently=True)
    logger.info('- email_notify_task_completion')
    return

class ValidateTaskView(View):
    """View to handle dynamic loading of validation results from `viewer.tasks.validate`.
    The validation of files uploaded to viewer/upload_cset or a target set by a user
    at viewer/upload_tset
    """
    def get(self, request, validate_task_id):
        """Get method for `ValidateTaskView`. Takes a validate task id, checks its
        status and returns the status, and result if the task is complete

        Returns
        -------
        response_data: JSON
            response data (dict) in JSON format:
                - if status = 'RUNNING':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                - if status = 'FAILURE':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                    - validate_traceback (str): task.traceback
                - if status = 'SUCCESS':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                    - html (str): html of task outcome - success message or html table of errors & fail message

        """
        logger.info('+ ValidateTaskView.get')
        validate_task_id_str = str(validate_task_id)

        task = AsyncResult(validate_task_id_str)
        response_data = {'validate_task_status': task.status,
                         'validate_task_id': task.id}

        if task.status == 'FAILURE':
            logger.info('+ ValidateTaskView.get.FAILURE')
            result = task.traceback
            response_data['validate_traceback'] = str(result)

            return JsonResponse(response_data)

        # Check if results ready
        if task.status == "SUCCESS":
            logger.info('+ ValidateTaskView.get.SUCCESS')
            results = task.get()
            # NB get tuple from validate task
            process_type = results[1]
            validate_dict = results[2]
            validated = results[3]
            if validated:
                response_data['html'] = 'Your data was validated. \n It can now be uploaded using the upload option.'
                response_data['validated'] = 'Validated'

                if process_type== 'tset':
                    target_name = results[5]
                    contact_email = results[8]
                    email_task_completion(contact_email, 'validate-success', target_name)

                return JsonResponse(response_data)

            if not validated:
                # set pandas options to display all column data
                pd.set_option('display.max_colwidth', -1)

                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''

                response_data["html"] = html_table
                response_data['validated'] = 'Not validated'
                if process_type== 'tset':
                    target_name = results[5]
                    contact_email = results[8]
                    email_task_completion(contact_email, 'validate-failure', target_name, task_id=validate_task_id_str)

                return JsonResponse(response_data)

        return JsonResponse(response_data)


class UpdateTaskView(View):

    def get(self, request, update_task_id):
        update_task_id_str = str(update_task_id)
        task = AsyncResult(update_task_id_str)
        response_data = {'update_task_status': task.status,
                         'update_task_id': task.id}

        result = 'Running...'

        if task.status == 'FAILURE':
            result = task.traceback
            response_data['result'] = str(result)

        if task.status == 'SUCCESS':
            result = task.get()

        response_data['result'] = str(result)

        return JsonResponse(response_data)


class UploadTaskView(View):
    """View to handle dynamic loading of upload results from `viewer.tasks.process_compound_set`.
    The upload of files for a computed set by a user at viewer/upload_cset or a target
    set by a user at viewer/upload_tset.
    """
    def get(self, request, upload_task_id):
        """Get method for `UploadTaskView`. Takes an upload task id, checks its
        status and returns the status, and result if the task is complete

        Returns
        -------
        response_data: JSON
            response data (dict) in JSON format:
                - if status = 'RUNNING':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                - if status = 'FAILURE':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                    - upload_traceback (str): task.traceback
                - if status = 'SUCCESS':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                    - if results are a list (data was processed - validated or uploaded):
                        if this was a validation process
                        - validated (str): 'Not validated'
                        - html (str): html table of validation errors
                        if results are a validation/upload process:
                        - validated (str): 'Validated'
                        - results (dict): results
                        For compound sets ('cset')
                        - results['cset_download_url'] (str): download url for computed set sdf file
                        - results['pset_download_url'] (str): download url for computed set pdb files (zip)
                        For target sets ('tset')
                        - results['tset_download_url'] (str): download url for processed zip file
                    - if results are not string or list:
                        - processed (str): 'None'
                        - html (str): message to tell the user their data was not processed
        """
        logger.debug('+ UploadTaskView.get')
        upload_task_id_str = str(upload_task_id)
        task = AsyncResult(upload_task_id_str)
        response_data = {'upload_task_status': task.status,
                         'upload_task_id': task.id}

        if task.status == 'FAILURE':
            result = task.traceback
            response_data['upload_traceback'] = str(result)

            return JsonResponse(response_data)

        if task.status == 'SUCCESS':
            logger.debug('+ UploadTaskView.get.success')

            results = task.get()

            # Validation output for a cset or tset is a dictionary.
            if isinstance(results, list):
                if results[0] == 'validate':
                    # Get dictionary results
                    validate_dict = results[1]

                    # set pandas options to display all column data
                    pd.set_option('display.max_colwidth', -1)
                    table = pd.DataFrame.from_dict(results[2])
                    html_table = table.to_html()
                    html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''

                    response_data['validated'] = 'Not validated'
                    response_data['html'] = html_table

                    return JsonResponse(response_data)
                else:
                    # Upload/Update output tasks send back a tuple
                    # First element defines the source of the upload task (cset, tset)
                    response_data['validated'] = 'Validated'
                    if results[1] == 'tset':
                        target_name = results[2]
                        contact_email = results[5]
                        target_path = '/viewer/target/%s' % target_name
                        response_data['results'] = {}
                        response_data['results']['tset_download_url'] = target_path
                        logger.info('+ UploadTaskView.get.success -email: %s', contact_email)
                        email_task_completion(contact_email, 'upload-success', target_name, target_path=target_path)
                    else:
                        cset_name = results[2]
                        cset = models.ComputedSet.objects.get(name=cset_name)
                        submitter = cset.submitter
                        name = cset.unique_name
                        response_data['results'] = {}
                        response_data['results']['cset_download_url'] = '/viewer/compound_set/%s' % name
                        response_data['results']['pset_download_url'] = '/viewer/protein_set/%s' % name

                    return JsonResponse(response_data)

            else:
                # Error output
                html_table = '''<p> Your data was <b>not</b> processed.</p>'''
                response_data['processed'] = 'None'
                response_data['html'] = html_table
                return JsonResponse(response_data)

        return JsonResponse(response_data)


def img_from_smiles(request):
    """Generate a 2D molecule image for a given smiles string
    """
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        if smiles:
            return get_params(smiles, request)
        else:
            return HttpResponse("Please insert SMILES")
    else:
        return HttpResponse("Please insert SMILES")


def highlight_mol_diff(request):
    """Generate a 2D molecule image highlighting the difference between a
    reference and new molecule
    """
    if 'prb_smiles' and 'ref_smiles' in request.GET:
        return HttpResponse(get_highlighted_diffs(request))
    else:
        return HttpResponse("Please insert smiles for reference and probe")


def similarity_search(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
    else:
        return HttpResponse("Please insert SMILES")
    if "db_name" in request.GET:
        db_name = request.GET["db_name"]
    else:
        return HttpResponse("Please insert db_name")
    sql_query = """SELECT sub.*
  FROM (
    SELECT rdk.id,rdk.structure,rdk.idnumber
      FROM vendordbs.enamine_real_dsi_molfps AS mfp
      JOIN vendordbs.enamine_real_dsi AS rdk ON mfp.id = rdk.id
      WHERE m @> qmol_from_smiles(%s) LIMIT 1000
  ) sub;"""
    with connections[db_name].cursor() as cursor:
        cursor.execute(sql_query, [smiles])
        return HttpResponse(json.dumps(cursor.fetchall()))


def get_open_targets(request):
    """Return a list of all open targets (viewer/open_targets)
    """
    targets = models.Target.objects.all()
    target_names = []
    target_ids = []

    for t in targets:
        for p in t.project_id.all():
            if 'OPEN' in p.title:
                target_names.append(t.title)
                target_ids.append(t.id)

    return HttpResponse(json.dumps({'target_names': target_names, 'target_ids': target_ids}))


# This is used in the URL on the process results page after uploading a compound_set
def cset_download(request, name):
    """View to download an SDF file of a computed set by name
    (viewer/compound_set/(<name>)).
    """
    compound_set = models.ComputedSet.objects.get(unique_name=name)
    filepath = compound_set.submitted_sdf
    with open(filepath.path, 'r', encoding='utf-8') as fp:
        data = fp.read()
    filename = 'compund-set_' + name + '.sdf'
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s' % filename  # force browser to download file
    response.write(data)
    return response


def pset_download(request, name):
    """View to download a zip file of all protein structures (apo) for a computed set
    (viewer/compound_set/(<name>))
     """
    response = HttpResponse(content_type='application/zip')
    filename = 'protein-set_' + name + '.zip'
    response['Content-Disposition'] = 'filename=%s' % filename  # force browser to download file

    compound_set = models.ComputedSet.objects.get(unique_name=name)
    computed = models.ComputedMolecule.objects.filter(computed_set=compound_set)
    pdb_filepaths = list(set([c.pdb_info.path for c in computed]))

    buff = StringIO()
    zip_obj = zipfile.ZipFile(buff, 'w')

    for fp in pdb_filepaths:
        data = open(fp, 'r', encoding='utf-8').read()
        zip_obj.writestr(fp.split('/')[-1], data)
    zip_obj.close()

    buff.flush()
    ret_zip = buff.getvalue()
    buff.close()
    response.write(ret_zip)

    return response


# This is used in the URL on the process results page after uploading a target_set
def tset_download(request, title):
    """View to download an zip file of a target set by name (viewer/target/(<title>)).
    """
    target_set = models.Target.objects.get(title=title)
    media_root = settings.MEDIA_ROOT
    filepath = os.path.join(media_root, target_set.zip_archive.name)
    target_zip = open(filepath, 'rb')
    filename = 'target-set_' + title + '.zip'
    response = HttpResponse(target_zip, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename="%s"' % filename  # force browser to download file
    return response


# Start of ActionType
class ActionTypeView(viewsets.ModelViewSet):
    """View to retrieve information about action types available to users (GET).
    (api/action-types).
    """
    queryset = models.ActionType.objects.filter()
    serializer_class = serializers.ActionTypeSerializer

    # POST method allowed for flexibility in the PoC. In the final design we may want to prevent the POST/PUT methods
    # from being used
    # for action types so that these can only be updated via the admin panel.
    #    http_method_names = ['get', 'head']

    filterset_fields = '__all__'


# Start of Session Project
class SessionProjectsView(viewsets.ModelViewSet):
    """View to retrieve information about user projects (collection of sessions) (GET).
    Also used for saving project information (PUT, POST, PATCH).
    (api/session-projects).
    """
    queryset = models.SessionProject.objects.filter()

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a POST, PUT or PATCH request

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.SessionProjectReadSerializer`
            - if other: `viewer.serializers.SessionProjectWriteSerializer`
        """
        if self.request.method in ['GET']:
            # GET
            return serializers.SessionProjectReadSerializer
        # (POST, PUT, PATCH)
        return serializers.SessionProjectWriteSerializer

    filter_permissions = "target_id__project_id"
    filterset_fields = '__all__'


class SessionActionsView(viewsets.ModelViewSet):
    """View to retrieve information about actions relating to sessions_project (GET).
    Also used for saving project action information (PUT, POST, PATCH).
    (api/session-actions).
     """
    queryset = models.SessionActions.objects.filter()
    serializer_class = serializers.SessionActionsSerializer

    #   Note: jsonField for Actions will need specific queries - can introduce if needed.
    filterset_fields = ('id', 'author', 'session_project', 'last_update_date')


class SnapshotsView(viewsets.ModelViewSet):
    """View to retrieve information about user sessions (snapshots) (GET).
    Also used for saving session information (PUT, POST, PATCH). (api/snapshots)
    """
    queryset = models.Snapshot.objects.filter()

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a POST, PUT or PATCH request

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.SnapshotReadSerializer`
            - if other: `viewer.serializers.SnapshotWriteSerializer`
        """
        if self.request.method in ['GET']:
            return serializers.SnapshotReadSerializer
        return serializers.SnapshotWriteSerializer

    filter_class = filters.SnapshotFilter


class SnapshotActionsView(viewsets.ModelViewSet):
    """View to retrieve information about actions relating to snapshots (GET).
    Also used for saving snapshot action information (PUT, POST, PATCH).
    (api/snapshot-actions).
    """
    queryset = models.SnapshotActions.objects.filter()
    serializer_class = serializers.SnapshotActionsSerializer

    #   Note: jsonField for Actions will need specific queries - can introduce if needed.
    filterset_fields = ('id', 'author', 'session_project', 'snapshot', 'last_update_date')


class DSetCSVParser(BaseParser):
    """
    CSV parser class specific to design set csv spec -
    sets media_type for DSetUploadView to text/csv
    """
    media_type = 'text/csv'


class DSetUploadView(APIView):
    """Upload a design set (PUT) from a csv file
    """
    parser_class = (DSetCSVParser,)

    def put(self, request, format=None):  # pylint: disable=redefined-builtin
        """Method to handle PUT request and upload a design set
        """
        # Don't need...
        del format

        f = request.FILES['file']
        set_type = request.PUT['type']
        set_description = request.PUT['description']

        # save uploaded file to temporary storage
        name = f.name
        path = default_storage.save('tmp/' + name, ContentFile(f.read()))
        tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

        df = pd.read_csv(tmp_file)
        mandatory_cols = ['set_name', 'smiles', 'identifier', 'inspirations']
        actual_cols = df.columns
        for col in mandatory_cols:
            if col not in actual_cols:
                raise ParseError("The 4 following columns are mandatory: set_name, smiles, identifier, inspirations")

        set_names, compounds = process_design_sets(df, set_type, set_description)

        string = 'Design set(s) successfully created: '

        length = len(set_names)
        string += str(length) + '; '
        for i in range(0, length):
            string += str(i + 1) + ' - ' + set_names[i] + ') number of compounds = ' + str(len(compounds[i])) + '; '

        return HttpResponse(json.dumps(string))


class ComputedSetView(viewsets.ModelViewSet):
    """Retrieve information about and delete computed sets.
    """
    queryset = models.ComputedSet.objects.filter()
    serializer_class = serializers.ComputedSetSerializer
    filter_permissions = "project_id"
    filterset_fields = ('target', 'target__title')

    http_method_names = ['get', 'head', 'delete']

    def destroy(self, request, pk=None):
        """User provides the name of the ComputedSet (that's its primary key).
        We simply look it up and delete it, returning a standard 204 on success.
        """
        computed_set = get_object_or_404(models.ComputedSet, pk=pk)
        computed_set.delete()
        return HttpResponse(status=204)


class ComputedMoleculesView(viewsets.ReadOnlyModelViewSet):
    """Retrieve information about computed molecules - 3D info (api/compound-molecules).
    """
    queryset = models.ComputedMolecule.objects.filter()
    serializer_class = serializers.ComputedMoleculeSerializer
    filter_permissions = "project_id"
    filterset_fields = ('computed_set',)


class NumericalScoresView(viewsets.ReadOnlyModelViewSet):
    """View to retrieve information about numerical computed molecule scores
    (api/numerical-scores).
    """
    queryset = models.NumericalScoreValues.objects.filter()
    serializer_class = serializers.NumericalScoreSerializer
    filter_permissions = "project_id"
    filterset_fields = ('compound', 'score')


class TextScoresView(viewsets.ReadOnlyModelViewSet):
    """View to retrieve information about text computed molecule scores (api/text-scores).
    """
    queryset = models.TextScoreValues.objects.filter()
    serializer_class = serializers.TextScoreSerializer
    filter_permissions = "project_id"
    filterset_fields = ('compound', 'score')


class CompoundScoresView(viewsets.ReadOnlyModelViewSet):
    """View to retrieve descriptions of scores for a given name or computed set.
    """
    queryset = models.ScoreDescription.objects.filter()
    serializer_class = serializers.ScoreDescriptionSerializer
    filter_permissions = "project_id"
    filterset_fields = ('computed_set', 'name')


class ComputedMolAndScoreView(viewsets.ReadOnlyModelViewSet):
    """View to retrieve all information about molecules from a computed set
    along with all of their scores.
    """
    queryset = models.ComputedMolecule.objects.filter()
    serializer_class = serializers.ComputedMolAndScoreSerializer
    filter_permissions = "project_id"
    filterset_fields = ('computed_set',)


class DiscoursePostView(viewsets.ViewSet):
    """View to get and post to the Discourse platform

    example of input (GET) on local:

        http://127.0.0.1:8080/api/discourse_post/?post_title=Mpro%20First%20Project

    examples of input (POST raw data):

        {"category_name": "NewCategory", "parent_category_name": "Fragalysis targets", "category_colour": "0088CC",
        "category_text_colour": "FFFFFF", "post_title": "", "post_content": "", "post_tags":""}

        {"category_name": "NewCategory", "parent_category_name": "Fragalysis targets", "category_colour": "0088CC",
        "category_text_colour": "FFFFFF", "post_title": "New Topic Title 1",
        "post_content": "This is the first post that creates the topic - must be greater than 20 chars",
        "post_tags" :"[\\"tag1\\",\\"tag2\\"]"}

        {"category_name": "NewCategory", "parent_category_name": "Fragalysis targets", "category_colour": "0088CC",
        "category_text_colour": "FFFFFF", "post_title": "New Topic Title 1",
        "post_content": "This is a second post to New Topic Title 1", "post_tags" :"[]"}

    """
    serializer_class = serializers.DiscoursePostWriteSerializer

    def create(self, request):
        """Method to handle POST request and call discourse to create the post
        """
        logger.info('+ DiscoursePostView.post')
        data = request.data

        logger.info('+ DiscoursePostView.post %s', json.dumps(data))
        if data['category_name'] == '':
            category_details = None
        else:
            category_details = {'category_name': data['category_name'],
                                'parent_name': data['parent_category_name'],
                                'category_colour': data['category_colour'],
                                'category_text_colour': data['category_text_colour']}

        if data['post_title'] == '':
            post_details = None
        else:
            post_details = {'title': data['post_title'],
                            'content': data['post_content'],
                            'tags': json.loads(data['post_tags'])}

        error, post_url, error_message = create_discourse_post(request.user, category_details, post_details)

        logger.info('- DiscoursePostView.post')
        if error:
            return Response({"message": error_message})
        else:
            return Response({"Post url": post_url})

    def list(self, request):
        """Method to handle GET request and call discourse to list posts for a topic
        """
        logger.info('+ DiscoursePostView.get')
        query_params = request.query_params
        logger.info('+ DiscoursePostView.get %s', json.dumps(query_params))

        discourse_api_key = settings.DISCOURSE_API_KEY

        if discourse_api_key:
            post_title = request.query_params.get('post_title', None)
            error, posts = list_discourse_posts_for_topic(post_title)
        else:
            logger.info('- DiscoursePostView.no key')
            return Response({"message": "Discourse Not Available - No API key supplied"})

        logger.info('- DiscoursePostView.get')
        if error:
            return Response({"message": "No Posts Found"})
        else:
            return Response({"Posts": posts})


def create_csv_from_dict(input_dict, title=None, filename=None):
    """Write a CSV file containing data from an input dictionary and return a URL
    to the file in the media directory.
    """
    if not filename:
        filename = 'download'

    media_root = settings.MEDIA_ROOT
    unique_dir = str(uuid.uuid4())
    # /code/media/downloads/unique_dir
    download_path = os.path.join(media_root, 'downloads', unique_dir)
    os.makedirs(download_path, exist_ok=True)

    download_file = os.path.join(download_path, filename)

    # Remove file if it already exists
    if os.path.isfile(download_file):
        os.remove(download_file)

    with open(download_file, "w", newline='', encoding='utf-8') as csvfile:
        if title:
            csvfile.write(title)
            csvfile.write("\n")

    df = pd.DataFrame.from_dict(input_dict)
    df.to_csv(download_file, mode='a', header=True, index=False)

    return download_file


class DictToCsv(viewsets.ViewSet):
    """Takes a dictionary and returns a download link to a CSV file with the data.
    """
    serializer_class = serializers.DictToCsvSerializer

    def list(self, request):
        """Method to handle GET request
        """
        file_url = request.GET.get('file_url')

        if file_url and os.path.isfile(file_url):
            with open(file_url, encoding='utf8') as csvfile:
                # return file and tidy up.
                response = HttpResponse(csvfile, content_type='text/csv')
                response['Content-Disposition'] = 'attachment; filename=download.csv'
                shutil.rmtree(os.path.dirname(file_url), ignore_errors=True)
                return response
        else:
            return Response("Please provide file_url parameter")

    def create(self, request):
        """Method to handle POST request
        """
        logger.info('+ DictToCsv.post')
        input_dict = request.data['dict']
        input_title = request.data['title']

        if not input_dict:
            return Response({"message": "Please enter Dictionary"})
        else:
            filename_url = create_csv_from_dict(input_dict, input_title)

        return Response({"file_url": filename_url})


# Classes Relating to Tags
class TagCategoryView(viewsets.ModelViewSet):
    """Set up and retrieve information about tag categories (api/tag_category).
    """
    queryset = models.TagCategory.objects.filter()
    serializer_class = serializers.TagCategorySerializer
    filterset_fields = ('id', 'category')


class SiteObservationTagView(viewsets.ModelViewSet):
    """Set up/retrieve information about tags relating to Molecules (api/molecule_tag)
    """
    queryset = models.SiteObservationTag.objects.filter()
    serializer_class = serializers.SiteObservationTagSerializer
    filterset_fields = ('id', 'tag', 'category', 'target', 'molecules', 'mol_group')


class SessionProjectTagView(viewsets.ModelViewSet):
    """Set up/retrieve information about tags relating to Session Projects.
    """
    queryset = models.SessionProjectTag.objects.filter()
    serializer_class = serializers.SessionProjectTagSerializer
    filterset_fields = ('id', 'tag', 'category', 'target', 'session_projects')


class TargetMoleculesView(ISpyBSafeQuerySet):
    """Retrieve all Molecules and Tag information relating
    to a Target. The idea is that a single call can return all target related
    information needed by the React front end in a single call.
    """
    queryset = models.Target.objects.filter()
    serializer_class = serializers.TargetMoleculesSerializer
    filter_permissions = "project_id"
    filterset_fields = ("title",)


class DownloadStructures(ISpyBSafeQuerySet):
    """Uses a selected subset of the target data
    (proteins and booleans with suggested files) and creates a Zip file
    with the contents.

    Note that old zip files are removed after one hour.
    """
    queryset = models.Target.objects.filter()
    serializer_class = serializers.DownloadStructuresSerializer
    filter_permissions = "project_id"
    filterset_fields = ('title','id')

    def list(self, request):
        """Method to handle GET request
        """
        file_url = request.GET.get('file_url')

        if file_url:
            link = models.DownloadLinks.objects.filter(file_url=file_url)
            if (link and link[0].zip_file
                    and os.path.isfile(link[0].file_url)):
                logger.info('zip_file: %s', link[0].zip_file)

                # return file and tidy up.
                file_name = os.path.basename(file_url)
                wrapper = FileWrapper(open(file_url, 'rb'))
                response = FileResponse(wrapper,
                                        content_type='application/zip')
                response[
                    'Content-Disposition'] = \
                    'attachment; filename="%s"' % file_name
                response['Content-Length'] = os.path.getsize(file_url)
                return response
            elif link:
                content = {'message': 'Zip file no longer present - '
                                      'please recreate by calling '
                                      'POST/Prepare download'}
                return Response(content, status=status.HTTP_404_NOT_FOUND)

            content = {'message': 'File_url is not found'}
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        content = {'message': 'Please provide file_url parameter from '
                              'post response'}
        return Response(content, status=status.HTTP_404_NOT_FOUND)

    def create(self, request):
        """Method to handle POST request that's use to initiate a target download.\

        The user is permitted to download Targets they have ascess to (whether
        authenticated or not), and this is hadled by the queryset logic later in
        this method.
        """
        logger.info('+ DownloadStructures.post')

        # Clear up old existing files
        maintain_download_links()

        # Static files
        # For static files, the contents of the zip file at the time of the search
        # are stored in the zip_contents field. These are used to reconstruct the
        # zip file from the time of the request.
        if request.data['file_url']:
            # This is a static link - the contents are stored in the database
            # if required.
            file_url = request.data['file_url']
            logger.info('Given file_url "%s"', file_url)
            existing_link = models.DownloadLinks.objects.filter(file_url=file_url)

            if existing_link and existing_link[0].static_link:
                # If the zip file is present, return it
                # Note that don't depend 100% on the zip_file flag as the
                # file might have been deleted from the media server.
                if (existing_link[0].zip_file and
                        os.path.isfile(existing_link[0].file_url)):
                    logger.info('Download is Ready!')
                    return Response({"file_url": existing_link[0].file_url},
                                    status=status.HTTP_200_OK)
                elif os.path.isfile(existing_link[0].file_url):
                    # If the file is there but zip_file is false, then it is
                    # probably being rebuilt by a parallel process.
                    logger.info('Download is under construction')
                    content = {'message': 'Zip being rebuilt - '
                                          'please try later'}
                    return Response(content,
                                    status=status.HTTP_208_ALREADY_REPORTED)
                else:
                    # Otherwise re-create the file.
                    logger.info('Recreating download...')
                    recreate_static_file (existing_link[0], request.get_host())
                    return Response({"file_url": existing_link[0].file_url},
                                    status=status.HTTP_200_OK)

            msg = 'file_url should only be provided for static files'
            logger.warning(msg)
            content = {'message': msg}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)

        # Dynamic files

        if 'target_name' not in request.data:
            content = {'message': 'If no file_url, a target_name (title) must be provided'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)

        target_name = request.data['target_name']
        target = None
        logger.info('Given target_name "%s"', target_name)

        # Check target_name is valid
        # (it should natch the title of an existing target)
        for targ in self.queryset:
            if targ.title == target_name:
                target = targ
                break

        if not target:
            msg = f'Either the Target "{target_name}" is not present or you are not permitted access it'
            logger.warning(msg)
            content = {'message': msg}
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        logger.info('Found Target record %r', target)
        site_obvs = models.SiteObservation.objects.none()
        proteins_list = []
        if request.data['proteins']:

            logger.info('Given Proteins in request')
            # Get first part of protein code
            proteins_list = [p.strip().split(":")[0]
                             for p in request.data['proteins'].split(',')]
            logger.info('Given %s Proteins %s', len(proteins_list), proteins_list)

            logger.info('Looking for SiteObservation records for given Proteins...')
            # Filter by protein codes
            for code_first_part in proteins_list:
                # prot = models.Protein.objects.filter(code__contains=code_first_part).values()
                # I don't see why I need to drop out of django objects here
                prot = models.SiteObservation.objects.filter(code__contains=code_first_part)
                if prot.exists():
                    site_obvs = prot[0:1]
                else:
                    logger.warning('Could not find SiteObservation record for "%s"', code_first_part)

        else:

            logger.info('Request had no Proteins')
            logger.info('Looking for Protein records for %r...', target)
            # proteins = models.Protein.objects.filter(target_id=target.id).values()
            site_obvs = models.SiteObservation.objects.filter(
                experiment__experiment_upload__target=target)

        if not site_obvs.exists():
            content = {'message': 'Please enter list of valid protein codes '
                                  'for target: {}, proteins: {}'
                .format(target.title, proteins_list) }
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        protein_repr = ""
        for protein in site_obvs:
            protein_repr += "%r " % protein
        logger.info('Collected %s Protein records: %r', site_obvs.count(), protein_repr)

        filename_url, file_exists = check_download_links(request,
                                                         target,
                                                         site_obvs)
        if file_exists:
            return Response({"file_url": filename_url})
        else:
            content = {'message': 'Zip being rebuilt - please try later'}
            return Response(content,
                            status=status.HTTP_208_ALREADY_REPORTED)


class UploadTargetExperiments(viewsets.ModelViewSet):
    serializer_class = serializers.TargetExperimentWriteSerializer
    permission_class = [permissions.IsAuthenticated]
    http_method_names = ('post',)

    def get_view_name(self):
        return "Upload Target Experiments"

    def create(self, request, *args, **kwargs):
        logger.info("+ UploadTargetExperiments.create called")
        del args, kwargs

        serializer = self.get_serializer_class()(data=request.data)
        if serializer.is_valid():
            # User must have access to the Project
            # and the Target must be in the Project.

            contact_email = serializer.validated_data['contact_email']
            proposal_ref = serializer.validated_data['proposal_ref']
            filename = serializer.validated_data['file']

            # I'm not creating ExperimentUpload object here, so I
            # suppose there's no point in keeping this as ModelViewSet

            tmpdir = Path(settings.MEDIA_ROOT).joinpath('tmp')
            tmpdir.mkdir(exist_ok=True)
            target_file = tmpdir.joinpath(filename.name)
            handle_uploaded_file(target_file, filename)

            task = task_load_target.delay(str(target_file), proposal_ref=proposal_ref,
                                          contact_email=contact_email, user_id=request.user.pk)
            logger.info("+ UploadTargetExperiments.create  got Celery id %s", task.task_id)

            url = reverse('viewer:task_status', kwargs={'task_id': task.task_id})
            # as it launches task, I think 202 is more appropriate
            return Response({'task_status_url': url}, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class TaskStatus(APIView):
    def get(self, request, task_id, *args, **kwargs):
        """Given a task_id (a string) we try to return the status of the task,
        trying to handle unknown tasks as best we can. Celery is happy to accept any
        string as a Task ID. To know it's a real task takes a lot of work, or you can
        simply interpret a 'PENDING' state as "unknown task".
        """
        logger.debug("+ TaskStatus.get called task_id='%s'", task_id)

        # task_id is (will be) a UUID, but Celery expects a string
        task_id_str = str(task_id)
        try:
            result = AsyncResult(task_id_str)
        except TimeoutError:
            error = {'error': 'Task result query timed out'}
            return Response(error, status=status.HTTP_408_REQUEST_TIMEOUT)

        if result.state == 'PENDING':
            error = {'error': 'Unknown task'}
            return Response(error, status=status.HTTP_400_BAD_REQUEST)

        # Extract messages (from task info)
        # Assuming the task has some info.
        messages = []
        if result.info:
            # messages = [k for k in result.info.get('description', [])]
            messages = result.info.get('description', [])

        data = {
            'task_id': result.id,
            'status': result.state,
            'ready': result.ready(),
            'successful': result.successful(),
            'failed': result.failed(),
            'messages': messages,
        }
        return JsonResponse(data)


class TargetExperimentUploads(viewsets.ModelViewSet):
    queryset = models.ExperimentUpload.objects.all()
    serializer_class = serializers.TargetExperimentReadSerializer
    permission_class = [permissions.IsAuthenticated]
    filterset_fields = ("target", "project")
    http_method_names = ('get',)


class SiteObservations(viewsets.ModelViewSet):
    queryset = models.SiteObservation.objects.all()
    serializer_class = serializers.SiteObservationReadSerializer
    permission_class = [permissions.IsAuthenticated]
    http_method_names = ('get',)


class CanonSites(viewsets.ModelViewSet):
    queryset = models.CanonSite.objects.all()
    serializer_class = serializers.CanonSiteReadSerializer
    permission_class = [permissions.IsAuthenticated]
    http_method_names = ('get',)


class CanonSiteConfs(viewsets.ModelViewSet):
    queryset = models.CanonSiteConf.objects.all()
    serializer_class = serializers.CanonSiteConfReadSerializer
    permission_class = [permissions.IsAuthenticated]
    http_method_names = ('get',)


class XtalformSites(viewsets.ModelViewSet):
    queryset = models.XtalformSite.objects.all()
    serializer_class = serializers.XtalformSiteReadSerializer
    permission_class = [permissions.IsAuthenticated]
    http_method_names = ('get',)


class JobFileTransferView(viewsets.ModelViewSet):
    """Squonk Job file transfer (api/job_file_transfer)
    """
    queryset = models.JobFileTransfer.objects.filter()
    filter_permissions = "target__project_id"
    filterset_fields = ('id', 'snapshot', 'target', 'user',
                        'squonk_project', 'transfer_status')

    def get_serializer_class(self):
        if self.request.method in ['GET']:
            return serializers.JobFileTransferReadSerializer
        # (POST, PUT, PATCH)
        return serializers.JobFileTransferWriteSerializer

    def create(self, request):
        """Method to handle POST request
        """
        logger.info('+ JobFileTransfer.post')
        # Only authenticated users can transfer files to sqonk
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can transfer files'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the Squonk2 agent is not configured
        sq2a_rv = _SQ2A.configured()
        if not sq2a_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Collect expected API parameters....
        access_id = request.data['access']  # The access ID/legacy Project record title
        target_id = request.data['target']
        snapshot_id = request.data['snapshot']
        session_project_id = request.data['session_project']

        logger.info('+ user="%s" (id=%s)', user.username, user.id)
        logger.info('+ access_id=%s', access_id)
        logger.info('+ target_id=%s', target_id)
        logger.info('+ snapshot_id=%s', snapshot_id)
        logger.info('+ session_project_id=%s', session_project_id)

        target = models.Target.objects.get(id=target_id)
        assert target

        # Check the user can use this Squonk2 facility.
        # To do this we need to setup a couple of API parameter objects.
        sq2a_common_params: CommonParams = CommonParams(user_id=user.id,
                                                        access_id=access_id,
                                                        session_id=session_project_id,
                                                        target_id=target_id)
        sq2a_send_params: SendParams = SendParams(common=sq2a_common_params,
                                                  snapshot_id=snapshot_id)
        sq2a_rv: Squonk2AgentRv = _SQ2A.can_send(sq2a_send_params)
        if not sq2a_rv.success:
            content = {f'You cannot do this ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Check the presense of the files expected to be transferred
        error, proteins, compounds = check_file_transfer(request)
        if error:
            return Response(error['message'], status=error['status'])

        # Create new file transfer job
        logger.info('+ Calling ensure_project() to get the Squonk2 Project...')

        # This requires a Squonk2 Project (created by the Squonk2Agent).
        # It may be an existing project, or it might be a new project.
        common_params = CommonParams(user_id=user.id,
                                     access_id=access_id,
                                     target_id=target_id,
                                     session_id=session_project_id)
        sq2_rv = _SQ2A.ensure_project(common_params)
        if not sq2_rv.success:
            msg = f'Failed to get/create a Squonk2 Project' \
                    f' for User "{user.username}", Access ID {access_id},' \
                    f' Target ID {target_id}, and SessionProject ID {session_project_id}.' \
                    f' Got "{sq2_rv.msg}".' \
                    ' Cannot continue'
            content = {'message': msg}
            logger.error(msg)
            return Response(content, status=status.HTTP_404_NOT_FOUND)

        # The Squonk2Project record is in the response msg
        squonk2_project_uuid = sq2_rv.msg.uuid
        squonk2_unit_name = sq2_rv.msg.unit.name
        squonk2_unit_uuid = sq2_rv.msg.unit.uuid
        logger.info('+ ensure_project() returned Project uuid=%s (unit="%s" unit_uuid=%s)',
                    squonk2_project_uuid, squonk2_unit_name, squonk2_unit_uuid)

        job_transfer = models.JobFileTransfer()
        job_transfer.user = request.user
        job_transfer.proteins = [p['code'] for p in proteins]
        job_transfer.compounds = [c['name'] for c in compounds]
        # We should use a foreign key,
        # but to avoid migration issues with the existing code
        # we continue to use the project UUID string field.
        job_transfer.squonk_project = squonk2_project_uuid
        job_transfer.target = models.Target.objects.get(id=target_id)
        job_transfer.snapshot = models.Snapshot.objects.get(id=snapshot_id)

        # The 'transfer target' (a sub-directory of the transfer root)
        # For example the root might be 'fragalysis-files'
        # and the target may be `CD44MMA` so the targets will be written to the
        # Squonk project at fragalysis-files/CD44MMA
        assert job_transfer.target
        assert job_transfer.target.title
        transfer_target = job_transfer.target.title
        logger.info('+ transfer_target=%s', transfer_target)

        job_transfer.transfer_status = 'PENDING'
        job_transfer.transfer_datetime = None
        job_transfer.transfer_progress = None
        job_transfer.save()

        # The root (in the Squonk project) where files will be written for this Job.
        # Something like "fragalysis-files/hjyx" for new transfers,
        # "fragalysis-files" for existing transfers
        if job_transfer.sub_path:
            transfer_root = os.path.join(settings.SQUONK2_MEDIA_DIRECTORY, job_transfer.sub_path)
        else:
            transfer_root = settings.SQUONK2_MEDIA_DIRECTORY
        logger.info('+ transfer_root=%s', transfer_root)

        logger.info('oidc_access_token')
        logger.info(request.session['oidc_access_token'])

        logger.info('+ Starting transfer (celery) (job_transfer.id=%s)...',
                    job_transfer.id)
        job_transfer_task = process_job_file_transfer.delay(request.session['oidc_access_token'],
                                                            job_transfer.id)

        content = {'id' : job_transfer.id,
                   'transfer_root': transfer_root,
                   'transfer_target': transfer_target,
                   'transfer_status': job_transfer.transfer_status,
                   'transfer_task_id': str(job_transfer_task)}
        return Response(content,
                        status=status.HTTP_200_OK)


class JobConfigView(viewsets.ReadOnlyModelViewSet):
    """Calls Squonk to get a requested job configuration
    """
    def list(self, request):
        """Method to handle GET request
        """
        query_params = request.query_params
        logger.info('+ JobConfigView.get: %s', json.dumps(query_params))

        # Only authenticated users can have squonk jobs
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can access squonk jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the squonk variables are not set!
        sqa_rv = _SQ2A.configured()
        if not sqa_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sqa_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        job_collection = request.query_params.get('job_collection', None)
        job_name = request.query_params.get('job_name', None)
        job_version = request.query_params.get('job_version', None)
        # User must provide collection, name and version
        if not job_collection or not job_name or not job_version:
            content = {'Please provide job_collection, job_name and job_version'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)

        content = get_squonk_job_config(request,
                                        job_collection=job_collection,
                                        job_name=job_name,
                                        job_version=job_version)

        return Response(content)


class JobOverrideView(viewsets.ModelViewSet):
    queryset = models.JobOverride.objects.all().order_by('-id')

    def get_serializer_class(self):
        if self.request.method in ['GET']:
            return serializers.JobOverrideReadSerializer
        return serializers.JobOverrideWriteSerializer

    def create(self, request):
        logger.info('+ JobOverride.post')
        # Only authenticated users can transfer files to sqonk
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can provide Job overrides'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Override is expected to be a JSON string,
        # but protect against format issues
        override = request.data['override']
        try:
            override_dict = json.loads(override)
        except ValueError:
            content = {'error': 'The override is not valid JSON'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)
        # We could use a schema but that's comlex.
        # For now let's just insist on some key fields in the provided override: -
        # - global
        # - fragalysis-jobs
        if "global" not in override_dict:
            content = {'error': 'The override does not contain a "global" key'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)
        if "fragalysis-jobs" not in override_dict:
            content = {'error': 'The override does not contain a "fragalysis-jobs" key'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)
        if type(override_dict["fragalysis-jobs"]) != list:
            content = {'error': 'The override "fragalysis-jobs" key is not a list'}
            return Response(content, status=status.HTTP_400_BAD_REQUEST)

        job_override = models.JobOverride()
        job_override.override = override_dict
        job_override.author = user
        job_override.save()

        return Response({"id": job_override.id})


class JobRequestView(APIView):

    def get(self, request):
        logger.info('+ JobRequest.get')

        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can access squonk jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the Squonk2 agent is not configured
        sq2a_rv = _SQ2A.configured()
        if not sq2a_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Iterate through each record, for JobRequests that are not 'finished'
        # we call into Squonk to get an update. We then return the (possibly) updated
        # records to the caller.

        results = []
        snapshot_id = request.query_params.get('snapshot', None)

        if snapshot_id:
            logger.info('+ JobRequest.get snapshot_id=%s', snapshot_id)
            job_requests = models.JobRequest.objects.filter(snapshot=int(snapshot_id))
        else:
            logger.info('+ JobRequest.get snapshot_id=(unset)')
            job_requests = models.JobRequest.objects.all()

        for jr in job_requests:
            if not jr.job_has_finished():
                logger.info('+ JobRequest.get (id=%s) has not finished (job_status=%s)',
                            jr.id, jr.job_status)

                # Job's not finished, an opportunity to call into Squonk
                # To get the current status. To do this we'll need
                # the 'callback context' we supplied when launching the Job.
                logger.info('+ JobRequest.get (id=%s, code=%s) getting update from Squonk...',
                            jr.id, jr.code)
                sq2a_rv = _SQ2A.get_instance_execution_status(jr.code)
                # If the job's now finished, update the record.
                # If the call was successful we'll get None (not finished),
                # 'LOST', 'SUCCESS' or 'FAILURE'
                if not sq2a_rv.success:
                    logger.warning('+ JobRequest.get (id=%s, code=%s) check failed (%s)',
                                    jr.id, jr.code, sq2a_rv.msg)
                elif sq2a_rv.success and sq2a_rv.msg:
                    logger.info('+ JobRequest.get (id=%s, code=%s) new status is (%s)',
                                jr.id, jr.code, sq2a_rv.msg)
                    transition_time = str(datetime.utcnow())
                    transition_time_utc = parse(transition_time).replace(tzinfo=pytz.UTC)
                    jr.job_status = sq2a_rv.msg
                    jr.job_status_datetime = transition_time_utc
                    jr.job_finish_datetime = transition_time_utc
                    jr.save()
                else:
                    logger.info('+ JobRequest.get (id=%s, code=%s) is (probably) still running',
                                jr.id, jr.code)

            serializer = serializers.JobRequestReadSerializer(jr)
            results.append(serializer.data)

        num_results = len(results)
        logger.info('+ JobRequest.get num_results=%s', num_results)

        # Simulate the original paged API response...
        content = {'count': num_results,
                   'next': None,
                   'previous': None,
                   'results': results}
        return Response(content, status=status.HTTP_200_OK)

    def post(self, request):
        logger.info('+ JobRequest.post')
        # Only authenticated users can create squonk job requests.
        user = self.request.user
        if not user.is_authenticated:
            content = {'Only authenticated users can run jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the Squonk2 agent is not configured
        sq2a_rv = _SQ2A.configured()
        if not sq2a_rv.success:
            content = {f'The Squonk2 Agent is not configured ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Collect expected API parameters....
        target_id = request.data['target']
        snapshot_id = request.data['snapshot']
        session_project_id = request.data['session_project']
        access_id = request.data['access']  # The access ID/legacy Project record

        logger.info('+ user="%s" (id=%s)', user.username, user.id)
        logger.info('+ access_id=%s', access_id)
        logger.info('+ target_id=%s', target_id)
        logger.info('+ snapshot_id=%s', snapshot_id)
        logger.info('+ session_project_id=%s', session_project_id)

        # Check the user can use this Squonk2 facility.
        # To do this we need to setup a couple of API parameter objects.
        # We don't (at this point) care about the Job spec or callback URL.
        sq2a_common_params: CommonParams = CommonParams(user_id=user.id,
                                                        access_id=access_id,
                                                        session_id=session_project_id,
                                                        target_id=target_id)
        sq2a_run_job_params: RunJobParams = RunJobParams(common=sq2a_common_params,
                                                         job_spec=None,
                                                         callback_url=None)
        sq2a_rv: Squonk2AgentRv = _SQ2A.can_run_job(sq2a_run_job_params)
        if not sq2a_rv.success:
            content = {f'You cannot do this ({sq2a_rv.msg})'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        try:
            job_id, squonk_url_ext = create_squonk_job(request)
        except ValueError as error:
            logger.info('Job Request failed: %s', error)
            content = {'error': str(error)}
            return Response(content,
                            status=status.HTTP_400_BAD_REQUEST)

        logger.info('SUCCESS (job_id=%s squonk_url_ext=%s)', job_id, squonk_url_ext)

        content = {'id': job_id, 'squonk_url_ext': squonk_url_ext}
        return Response(content,
                        status=status.HTTP_200_OK)


class JobCallBackView(viewsets.ModelViewSet):
    """View to allow the Squonk system to update the status and job information for a
    specific job identified by a UUID.
    """
    queryset = models.JobRequest.objects.all()
    lookup_field = "code"
    http_method_names = ['get', 'head', 'put']

    def get_serializer_class(self):
        """Determine which serializer to use based on whether the request is a GET or a PUT

        Returns
        -------
        Serializer (rest_framework.serializers.ModelSerializer):
            - if GET: `viewer.serializers.JobCallBackWriteSerializer`
            - if other: `viewer.serializers.JobCallBackWriteSerializer`
        """
        if self.request.method in ['GET']:
            # GET
            return serializers.JobCallBackReadSerializer
        # PUT
        return serializers.JobCallBackWriteSerializer

    def update(self, request, code=None):
        """Response to a PUT on the Job-Callback.
        We're given a 'code' which we use to lookup the corresponding JobRequest
        (there'll only be one).
        """

        jr = models.JobRequest.objects.get(code=code)
        logger.info('+ JobCallBackView.update(code=%s) jr=%s', code, jr)

        # request.data is rendered as a dictionary
        if not request.data:
            return HttpResponse(status=204)

        j_status = request.data['job_status']
        # Get the appropriate SQUONK_STATUS...
        status_changed = False
        for squonk_status in models.JobRequest.SQUONK_STATUS:
            if squonk_status[0] == j_status and jr.job_status != j_status:
                jr.job_status = squonk_status[1]
                status_changed = True
                break

        if not status_changed:
            logger.info('+ JobCallBackView.update(code=%s) status=%s ignoring (no status change)', code, status)
            return HttpResponse(status=204)

        # This is now a chance to safely set the squonk_url_ext using the instance ID
        # present in the callback (if it's not already set). The instance is
        # placed at the end of the string, and is expected to be found
        # by process_compound_set_file() which loads the output file back
        # into Fragalysis
        if not jr.squonk_url_ext:
            jr.squonk_url_ext = create_squonk_job_request_url(request.data['instance_id'])
            logger.info("+ JobCallBackView.update(code=%s) jr.squonk_url_ext='%s'", code, jr.squonk_url_ext)

        # Update the state transition time,
        # assuming UTC.
        transition_time = request.data.get('state_transition_time')
        if not transition_time:
            transition_time = str(datetime.utcnow())
            logger.warning("+ JobCallBackView.update(code=%s) callback is missing state_transition_time"
                           " (using '%s')", code, transition_time)
        transition_time_utc = parse(transition_time).replace(tzinfo=pytz.UTC)
        jr.job_status_datetime = transition_time_utc

        logger.info('+ JobCallBackView.update(code=%s) status=%s transition_time=%s (new status)',
                    code, status, transition_time)

        # If the Job's start-time is not set, set it.
        if not jr.job_start_datetime:
            logger.info('+ JobCallBackView.update(code=%s) setting job START datetime (%s)',
                        code, transition_time)
            jr.job_start_datetime = transition_time_utc

        # Set the Job's finish time (once) if it looks lie the Job's finished.
        # We can assume the Job's finished if the status is one of a number
        # of values...
        if not jr.job_finish_datetime and status in ('SUCCESS', 'FAILURE', 'REVOKED'):
            logger.info('+ JobCallBackView.update(code=%s) Setting job FINISH datetime (%s)',
                        code, transition_time)
            jr.job_finish_datetime = transition_time_utc

        # Save the JobRequest record before going further.
        jr.save()

        if j_status != 'SUCCESS':
            # Go no further unless SUCCESS
            return HttpResponse(status=204)

        logger.info('+ JobCallBackView.update(code=%s) job finished (SUCCESS).'
                    ' Can we upload the results?', code)

        # SUCCESS ... automatic upload?
        #
        # Only continue if the target file is 'merged.sdf'.
        # For now there must be an '--outfile' in the job info's 'command'.
        # Here we have hard-coded the expectations because the logic to identify the
        # command's outputs is not fully understood.
        # The command is a string that we split and search.
        job_output = ''
        jr_job_info_msg = jr.squonk_job_info[1]
        command = jr_job_info_msg.get('command')
        command_parts = shlex.split(command)
        outfile_index = 0
        while outfile_index < len(command_parts)\
                and command_parts[outfile_index] != '--outfile':
            outfile_index += 1
        # Found '--command'?
        if command_parts[outfile_index] == '--outfile'\
                and outfile_index < len(command_parts) - 1:
            # Yes ... the filename is the next item in the list
            job_output = command_parts[outfile_index + 1]
        job_output_path = '/' + os.path.dirname(job_output)
        job_output_filename = os.path.basename(job_output)

        logging.info('+ JobCallBackView.update(code=%s) job_output_path="%s"', code, job_output_path)
        logging.info('+ JobCallBackView.update(code=%s) job_output_filename="%s"', code, job_output_filename)

        # If it's not suitably named, leave
        expected_squonk_filename = 'merged.sdf'
        if job_output_filename != expected_squonk_filename:
            # Incorrectly named file - nothing to get/upload.
            logger.info('+ JobCallBackView.update(code=%s) SUCCESS but not uploading.'
                        ' Expected "%s" as job_output_filename.'
                        ' Found "%s"', code, expected_squonk_filename, job_output_filename)
            return HttpResponse(status=204)

        if jr.upload_status != 'PENDING':
            logger.warning('+ JobCallBackView.update(code=%s) SUCCESS but ignoring.'
                           ' upload_status=%s (already uploading?)', code, jr.upload_status)
            return HttpResponse(status=204)

        # Change of status and SUCCESS
        # - mark the job upload as 'started'
        jr.upload_status = "STARTED"
        jr.save()

        # Initiate an upload (and removal) of files from Squonk.
        # Which requires the linking of several tasks.
        # We star the process with 'process_compound_set_job_file'
        # with the path and filename already discoverd...
        task_params = {'jr_id': jr.id,
                       'transition_time': transition_time,
                       'job_output_path': job_output_path,
                       'job_output_filename': job_output_filename}
        task_upload = (
             process_compound_set_job_file.s(task_params) |
             validate_compound_set.s() |
             process_compound_set.s() |
             erase_compound_set_job_material.s(job_request_id=jr.id)
        ).apply_async()

        logger.info('+ JobCallBackView.update(code=%s)'
                    ' started process_job_file_upload(%s) task_upload=%s',
                    code, jr.id, task_upload)

        return HttpResponse(status=204)

class JobAccessView(APIView):
    """JobAccess (api/job_access)

    Django view that calls Squonk to allow a user (who is able to see a Job)
    the ability to access that Job in Squonk. To be successful the user
    must have access to the corresponding Fragalysis Project. This can be called by
    the Job 'owner', who always has access.
    """
    def get(self, request):
        """Method to handle GET request
        """
        query_params = request.query_params
        logger.info('+ JobAccessView/GET %s', json.dumps(query_params))

        err_response = {'accessible': False}
        ok_response = {'accessible': True, 'error': ''}

        # Only authenticated users can have squonk jobs
        user = self.request.user
        if not user.is_authenticated:
            content = {'accessible': False,
                       'error': 'Only authenticated users can access Squonk Jobs'}
            return Response(content, status=status.HTTP_403_FORBIDDEN)

        # Can't use this method if the squonk variables are not set!
        sqa_rv = _SQ2A.configured()
        if not sqa_rv.success:
            err_response['error'] = f'Squonk is not available ({sqa_rv.msg})'
            return Response(err_response, status=status.HTTP_403_FORBIDDEN)

        # Get the JobRequest Record form the supplied record ID
        jr_id_str = request.query_params.get('job_request_id', None)
        if not jr_id_str or not jr_id_str.isdigit():
            err_response['error'] = f'The JobRequest ID ({jr_id_str}) is not valid'
            return Response(err_response, status=status.HTTP_400_BAD_REQUEST)
        jr_id = int(jr_id_str)
        if jr_id < 1:
            err_response['error'] = f'The JobRequest ID ({jr_id}) cannot be less than 1'
            return Response(err_response, status=status.HTTP_400_BAD_REQUEST)

        jr_list = models.JobRequest.objects.filter(id=jr_id)
        if len(jr_list) == 0:
            err_response['error'] = 'The JobRequest does not exist'
            return Response(err_response, status=status.HTTP_400_BAD_REQUEST)
        jr = jr_list[0]

        # JobRequest must have a Squonk Project value
        if not jr.squonk_project:
            err_response['error'] = f'The JobRequest ({jr_id}) has no Squonk Project value'
            return Response(err_response, status=status.HTTP_403_FORBIDDEN)

        # User must have access to the Job's Project.
        # If the user is not the owner of the Job, and there is a Project,
        # we then check the user has access to the given access ID.
        #
        # If the Job has no Project (Jobs created before this change will not have a Project)
        # or the user is the owner of the Job we skip this check.
        if user.id != jr.user.id:
            logger.info('+ JobAccessView/GET Checking access to JobRequest %s for "%s" (%s)',
                        jr_id, user.username, jr.squonk_project)
            if not jr.project or not jr.project.title:
                logger.warning('+ JobAccessView/GET No Fragalysis Project (or title)'
                               ' for JobRequest %s - granting access',
                               jr_id)
            else:
                # The project is the Job's access ID.
                # To check access we need this and the User's ID
                access_id = jr.project.id
                access_string = jr.project.title
                sq2a_common_params: CommonParams = CommonParams(user_id=user.id,
                                                                access_id=access_id,
                                                                session_id=None,
                                                                target_id=None)
                sq2a_run_job_params: RunJobParams = RunJobParams(common=sq2a_common_params,
                                                                 job_spec=None,
                                                                 callback_url=None)
                sq2a_rv: Squonk2AgentRv = _SQ2A.can_run_job(sq2a_run_job_params)
                if not sq2a_rv.success:
                    err_response['error'] = f'Access to the Job for {access_id} ({access_string}) is denied. {sq2a_rv.msg}'
                    return Response(err_response, status=status.HTTP_403_FORBIDDEN)

            # All checks passed ... grant access for this user
            sq2a_access_params: AccessParams = AccessParams(username=user.username,
                                                            project_uuid=jr.squonk_project)
            sqa_rv = _SQ2A.grant_access(sq2a_access_params)
            if not sqa_rv.success:
                err_response['error'] = f'Squonk failed to grant access ({sqa_rv.msg})'
                logger.warning('+ JobAccessView/GET error=%s', content['error'])
                return Response(err_response, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        logger.info('+ JobAccessView/GET Success for %s/"%s" on %s',
                    jr_id, user.username, jr.squonk_project)
        return Response(ok_response)
