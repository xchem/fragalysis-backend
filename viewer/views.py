import json
import os
import zipfile
from cStringIO import StringIO

from django.db import connections
from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.conf import settings
from django.http import JsonResponse

from rest_framework.parsers import JSONParser, BaseParser
from rest_framework.exceptions import ParseError
from rest_framework.response import Response
from rest_framework.views import APIView

from django.views import View

from celery import current_app, chain
from celery.result import AsyncResult

from api.security import ISpyBSafeQuerySet
from api.utils import get_params, get_highlighted_diffs

from viewer.models import (
    Molecule,
    Protein,
    Compound,
    Target,
    SessionProject,
    Snapshot,
    ComputedMolecule,
    ComputedSet,
    CSetKeys,
    NumericalScoreValues,
    ScoreDescription,
    File
)
from viewer import filters
from forms import CSetForm, UploadKeyForm

from tasks import *

import pandas as pd

from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
    MolImageSerialzier,
    CmpdImageSerialzier,
    ProtMapInfoSerialzer,
    ProtPDBInfoSerialzer,
    ProtPDBBoundInfoSerialzer,
    VectorsSerializer,
    GraphSerializer,
    SessionProjectWriteSerializer,
    SessionProjectReadSerializer,
    SnapshotReadSerializer,
    SnapshotWriteSerializer,
    FileSerializer,
    ComputedSetSerializer,
    ComputedMoleculeSerializer,
    NumericalScoreSerializer,
    ScoreDescriptionSerializer,
    TextScoreSerializer,
    ComputedMolAndScoreSerializer,
)


class VectorsView(ISpyBSafeQuerySet):
    queryset = Molecule.objects.filter()
    serializer_class = VectorsSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class GraphView(ISpyBSafeQuerySet):
    queryset = Molecule.objects.filter()
    serializer_class = GraphSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class MolImageView(ISpyBSafeQuerySet):
    queryset = Molecule.objects.filter()
    serializer_class = MolImageSerialzier
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class CompoundImageView(ISpyBSafeQuerySet):
    queryset = Compound.objects.filter()
    serializer_class = CmpdImageSerialzier
    filter_permissions = "project_id"
    filter_fields = ("smiles",)


class ProteinMapInfoView(ISpyBSafeQuerySet):
    queryset = Protein.objects.filter()
    serializer_class = ProtMapInfoSerialzer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


class ProteinPDBInfoView(ISpyBSafeQuerySet):
    queryset = Protein.objects.filter()
    serializer_class = ProtPDBInfoSerialzer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


class ProteinPDBBoundInfoView(ISpyBSafeQuerySet):
    queryset = Protein.objects.filter()
    serializer_class = ProtPDBBoundInfoSerialzer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


class TargetView(ISpyBSafeQuerySet):
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    filter_permissions = "project_id"
    filter_fields = ("title",)


class MoleculeView(ISpyBSafeQuerySet):
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = (
        "prot_id",
        "cmpd_id",
        "smiles",
        "prot_id__target_id",
        "mol_type",
        "mol_groups",
    )


class CompoundView(ISpyBSafeQuerySet):
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_permissions = "project_id"
    filter_fields = ("smiles",)


class ProteinView(ISpyBSafeQuerySet):
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


def react(request):
    """
    :param request:
    :return:
    """
    return render(request, "viewer/react_temp.html")


# email cset upload key
def cset_key(request):
    form = UploadKeyForm()
    if request.method == 'POST':
        form = UploadKeyForm(request.POST)
        email = request.POST['contact_email']
        new_key = CSetKeys()
        new_key.user = email
        new_key.save()
        key_value = new_key.uuid

        from django.core.mail import send_mail
        from django.conf import settings

        subject = 'Fragalysis: upload compound set key'
        message = 'Your upload key is: ' + str(
            key_value) + ' store it somewhere safe. Only one key will be issued per user'
        email_from = settings.EMAIL_HOST_USER
        recipient_list = [email, ]
        send_mail(subject, message, email_from, recipient_list)

        msg = 'Your key will be emailed to: <b>' + email + '</b>'

        return render(request, 'viewer/generate-key.html', {'form': form, 'message': msg})
    return render(request, 'viewer/generate-key.html', {'form': form, 'message': ''})


class UploadCSet(View):

    def get(self, request):

        #test = TargetView().get_queryset(request=request)
        #targets = request.get('/api/targets/')
        #int(targets)
        form = CSetForm()
        return render(request, 'viewer/upload-cset.html', {'form': form})

    def post(self, request):
        check_services()
        zfile = None
        zf = None
        cset = None
        form = CSetForm(request.POST, request.FILES)
        context = {}
        if form.is_valid():
            # get the upload key
            # key = request.POST['upload_key']
            # all_keys = CSetKeys.objects.all()
            # if it's not valid, return a message
            # if key not in [str(key.uuid) for key in all_keys]:
            #     html = "<br><p>You either didn't provide an upload key, or it wasn't valid. Please try again (email rachael.skyner@diamond.ac.uk to obtain an upload key)</p>"
            #     return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html})

            # get all of the variables needed from the form
            myfile = request.FILES['sdf_file']
            target = request.POST['target_name']
            choice = request.POST['submit_choice']

            if 'pdb_zip' in request.FILES.keys():
                pdb_file = request.FILES['pdb_zip']
            else:
                pdb_file = None

            # if there is a zip file of pdbs, check it for .pdb files, and ignore others
            if pdb_file:
                zf = zipfile.ZipFile(pdb_file)
                zip_lst = zf.namelist()
                zip_names = []
                for filename in zip_lst:
                    # only handle pdb files
                    if filename.split('.')[-1] == 'pdb':
                        # store filenames?
                        zip_names.append(filename)

                zfile = {'zip_obj': zf, 'zf_list': zip_names}

            # Close the zip file
            if zf:
                zf.close()

            # save uploaded sdf to tmp storage
            name = myfile.name
            path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
            tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

            # Settings for if validate option selected
            if str(choice) == '0':
                # Start celery task
                task_validate = validate.delay(tmp_file, target=target, zfile=zfile)

                context = {}
                context['validate_task_id'] = task_validate.id
                context['validate_task_status'] = task_validate.status

                # Update client side with task id and status
                return render(request, 'viewer/upload-cset.html', context)

            # if it's an upload, run the compound set task
            if str(choice) == '1':
                # Start chained celery tasks. NB first function passes tuple
                # to second function - see tasks.py
                task_upload = (validate.s(tmp_file, target=target, zfile=zfile) | process_compound_set.s()).apply_async()

                context = {}
                context['upload_task_id'] = task_upload.id
                context['upload_task_status'] = task_upload.status

                # Update client side with task id and status
                return render(request, 'viewer/upload-cset.html', context)

        context['form'] = form
        return render(request, 'viewer/upload-cset.html', context)


# Add ValidateTaskView
class ValidateTaskView(View):

    def get(self, request, validate_task_id):
        task = AsyncResult(validate_task_id)
        response_data = {'validate_task_status': task.status,
                         'validate_task_id': task.id}

        if task.status == 'FAILURE':
            result = task.traceback
            response_data['validate_traceback'] = str(result)

            return JsonResponse(response_data)

        # Check if results ready
        if task.status == "SUCCESS":
            results = task.get()
            # NB get tuple from validate task
            validate_dict = results[0]
            validated = results[1]
            if validated:
                response_data['html'] = 'Your data was validated. \n It can now be uploaded using the upload option.'

                return JsonResponse(response_data)


            if not validated:
                # set pandas options to display all column data
                pd.set_option('display.max_colwidth', -1)

                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''

                response_data["html"] = html_table

                return JsonResponse(response_data)

        return JsonResponse(response_data)


# Needs to be something like UploadTaskView
class UploadTaskView(View):
    def get(self, request, upload_task_id):
        task = AsyncResult(upload_task_id)
        response_data = {'upload_task_status': task.status,
                         'upload_task_id': task.id}

        if task.status == 'FAILURE':
            result = task.traceback
            response_data['upload_traceback'] = str(result)

            return JsonResponse(response_data)

        if task.status == 'SUCCESS':

            results = task.get()

            # Check for d,v vs csetname output
            if isinstance(results, list):
                # Get dictionary results
                validate_dict = results[0]

                # set pandas options to display all column data
                pd.set_option('display.max_colwidth', -1)

                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''

                response_data['validated'] = 'Not validated'
                response_data['html'] = html_table

                return JsonResponse(response_data)

            # Check for d,v vs csetname output
            # Check in with Rachael if we are expecting a string here?
            if isinstance(results, str):
                cset_name = results
                cset = ComputedSet.objects.get(name=cset_name)

                submitter = cset.submitter
                name = cset.unique_name
                response_data['validated'] = 'Validated'
                response_data['results'] = {}
                response_data['results']['cset_download_url'] = '/viewer/compound_set/%s' % name
                response_data['results']['pset_download_url'] = '/viewer/protein_set/%s' % name

                return JsonResponse(response_data)

            else:

                html_table = '''<p> Your data was <b>not</b> processed.</p>'''
                response_data['processed'] = 'None'
                response_data['html'] = html_table
                return JsonResponse(response_data)

        return JsonResponse(response_data)


# Add SaveFilesTaskView

def img_from_smiles(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        if smiles:
            return get_params(smiles, request)
        else:
            return HttpResponse("Please insert SMILES")
    else:
        return HttpResponse("Please insert SMILES")


def highlight_mol_diff(request):
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
    targets = Target.objects.all()
    target_names = []
    target_ids = []

    for t in targets:
        for p in t.project_id.all():
            if 'OPEN' in p.title:
                target_names.append(t.title)
                target_ids.append(t.id)

    return HttpResponse(json.dumps({'target_names': target_names, 'target_ids': target_ids}))


def cset_download(request, name):
    compound_set = ComputedSet.objects.get(unique_name=name)
    filepath = compound_set.submitted_sdf
    with open(filepath.path, 'r') as fp:
        data = fp.read()
    filename = 'compund-set_' + name + '.sdf'
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s' % filename  # force browser to download file
    response.write(data)
    return response


def pset_download(request, name):
    response = HttpResponse(content_type='application/zip')
    filename = 'protein-set_' + name + '.zip'
    response['Content-Disposition'] = 'filename=%s' % filename  # force browser to download file

    compound_set = ComputedSet.objects.get(unique_name=name)
    computed = ComputedMolecule.objects.filter(computed_set=compound_set)
    pdb_filepaths = list(set([c.pdb_info.path for c in computed]))

    buff = StringIO()
    zip_obj = zipfile.ZipFile(buff, 'w')

    for fp in pdb_filepaths:
        data = open(fp, 'r').read()
        zip_obj.writestr(fp.split('/')[-1], data)
    zip_obj.close()

    buff.flush()
    ret_zip = buff.getvalue()
    buff.close()
    response.write(ret_zip)

    return response


## Start of Session Project
class SessionProjectsView(viewsets.ModelViewSet):
    queryset = SessionProject.objects.filter()

    def get_serializer_class(self):
        if self.request.method in ['GET']:
            # GET
            return SessionProjectReadSerializer
        # (POST, PUT, PATCH)
        return SessionProjectWriteSerializer

    filter_permissions = "target_id__project_id"
    filter_fields = '__all__'


class SnapshotsView(viewsets.ModelViewSet):
    queryset = Snapshot.objects.filter()

    def get_serializer_class(self):
        if self.request.method in ['GET']:
            return SnapshotReadSerializer
        return SnapshotWriteSerializer

    filter_class = filters.SnapshotFilter

### End of Session Project


### Design sets upload
# Custom parser class for a csv file
class DSetCSVParser(BaseParser):
    """
    CSV parser class specific to design set csv spec
    """
    media_type = 'text/csv'


class DSetUploadView(APIView):
    parser_class = (DSetCSVParser,)

    def put(self, request, format=None):

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


class ComputedSetView(viewsets.ReadOnlyModelViewSet):
    queryset = ComputedSet.objects.filter()
    serializer_class = ComputedSetSerializer
    filter_permissions = "project_id"
    filter_fields = ('target',)

class ComputedMoleculesView(viewsets.ReadOnlyModelViewSet):
    queryset = ComputedMolecule.objects.filter()
    serializer_class = ComputedMoleculeSerializer
    filter_permissions = "project_id"
    filter_fields = ('computed_set',)


class NumericalScoresView(viewsets.ReadOnlyModelViewSet):
    queryset = NumericalScoreValues.objects.filter()
    serializer_class = NumericalScoreSerializer
    filter_permissions = "project_id"
    filter_fields = ('compound', 'score')


class TextScoresView(viewsets.ReadOnlyModelViewSet):
    queryset = TextScoreValues.objects.filter()
    serializer_class = TextScoreSerializer
    filter_permissions = "project_id"
    filter_fields = ('compound', 'score')


class CompoundScoresView(viewsets.ReadOnlyModelViewSet):
    queryset = ScoreDescription.objects.filter()
    serializer_class = ScoreDescriptionSerializer
    filter_permissions = "project_id"
    filter_fields = ('computed_set',)


class ComputedMolAndScoreView(viewsets.ReadOnlyModelViewSet):
    queryset = ComputedMolecule.objects.filter()
    serializer_class = ComputedMolAndScoreSerializer
    filter_permissions = "project_id"
    filter_fields = ('computed_set',)

