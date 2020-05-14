import json, os
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

from django.views import View

from celery import current_app

from api.security import ISpyBSafeQuerySet
from api.utils import get_params, get_highlighted_diffs

from viewer.models import Molecule, Protein, Compound, Target, SessionProject, Snapshot, ComputedCompound, CompoundSet, CSetKeys
from viewer import filters
from forms import CSetForm, UploadKeyForm
from tasks import check_services, process_compound_set, validate
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
    SnapshotWriteSerializer
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
        message = 'Your upload key is: ' + str(key_value) + ' store it somewhere safe. Only one key will be issued per user'
        email_from = settings.EMAIL_HOST_USER
        recipient_list = [email, ]
        send_mail(subject, message, email_from, recipient_list)

        msg = 'Your key will be emailed to: <b>' + email + '</b>'

        return render(request, 'viewer/generate-key.html', {'form': form, 'message':msg})
    return render(request, 'viewer/generate-key.html', {'form': form, 'message': ''})


# overall view for upload compound set - needs to include task for:
# - saving the data from pdb zip to temporary storage
# - validation
# - upload
# Also need to check that all data is going into the database
# Worth looking at chaining tasks together: https://docs.celeryproject.org/en/stable/userguide/canvas.html
# https://stackoverflow.com/questions/39099267/how-can-i-create-a-chain-of-conditional-subtasks-in-celery
class UploadCSet(View):
    def get(self, request):
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

            # where does this actually belong?
            if zf:
                zf.close()

            # save uploaded sdf to tmp storage
            name = myfile.name
            path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
            tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

            if str(choice) == '0':
                # validate the file (make this a task)
                # How to use two outputs in celery?
                #d, v = validate(tmp_file, target=target, zfile=zfile)

                task_validate = validate.delay(tmp_file, target=target, zfile=zfile)

                #d,v = task_validate.get()

                # Task id (get id)
                context['validate_task_id'] = d.id
                context['validate_task_status'] = d.status

                if v:
                    # need to add JS for handling validation task to template
                    return render(request, 'viewer/upload-cset.html', context)

                # if the data isn't validated make a table of errors and return it
                # All of this needs moving to the validate task view

                if not v:
                    # set pandas options to display all column data
                    pd.set_option('display.max_colwidth', -1)

                    table = pd.DataFrame.from_dict(d)
                    html_table = table.to_html()
                    html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''
                    return render(request, 'viewer/upload-cset.html',
                                  {'form': form, 'table': html_table, 'download_url': ''})
                    # return ValidationError('We could not validate this file')

            # if it's an upload, run the compound set task
            if str(choice) == '1':
                task_validate = validate.delay(tmp_file, target=target, zfile=zfile)

                d,v = task_validate.get()

                # If validation fails
                if not v:
                    # set pandas options to display all column data
                    pd.set_option('display.max_colwidth', -1)

                    table = pd.DataFrame.from_dict(d)
                    html_table = table.to_html()
                    html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''
                    return render(request, 'viewer/upload-cset.html',
                                  {'form': form, 'table': html_table, 'download_url': ''})

                if v:
                    # Move to class below
                    task_upload = process_compound_set.delay(target=target, filename=tmp_file, zfile=zfile)

                    cset_name = task_upload.get()
                    #cset_name = process_compound_set.delay(target=target, filename=tmp_file, zfile=zfile)

                    # Need to figure these out - same as above
                    #context['upload_task_id'] = cset_name.id
                    #context['upload_task_status'] = cset_name.status

                    return render(request, 'viewer/upload-cset.html', context)

        context['form'] = form
        return render(request, 'viewer/upload-cset.html', context)

# Add ValidateTaskView
class ValidateTaskView(View):
    def get(self, request, task_id):
        task = current_app.AsyncResult(task_id)
        response_data = {'validate_task_status': task.status, 'validate_task_id': task.id}

        if task.status == 'SUCCESS':
            validate_dict, validated = task.get()
            if validated:
                # return that the set was validated, and kick off upload task if requested
                pass
            else:
                # return the dict to be turned into table, and don't do upload
                pass

        return JsonResponse(response_data)

# Needs to be something like UploadTaskView
class UploadTaskView(View):
    def get(self, request, task_id):
        task = current_app.AsyncResult(task_id)
        response_data = {'upload_task_status': task.status, 'upload_task_id': task.id}

        if task.status == 'SUCCESS':
            cset_name = task.get()
            if cset_name:
                cset = CompoundSet.objects.get(name=cset_name)
                submitter = cset.submitter
                name = submitter.unique_name
                response_data['results'] = {}
                response_data['results']['cset_download_url'] = '/viewer/compound_set/%s' % name
                response_data['results']['pset_download_url'] = '/viewer/protein_set/%s' % name

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
    compound_set = CompoundSet.objects.get(submitter__unique_name=name)
    filepath = compound_set.submitted_sdf
    with open(filepath.path, 'r') as fp:
        data = fp.read()
    filename = 'compund-set_' + name + '.sdf'
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename=%s' % filename # force browser to download file
    response.write(data)
    return response


def pset_download(request, name):
    response = HttpResponse(content_type='application/zip')
    filename = 'protein-set_' + name + '.zip'
    response['Content-Disposition'] = 'filename=%s' % filename  # force browser to download file

    compound_set = CompoundSet.objects.get(submitter__unique_name=name)
    computed = ComputedCompound.objects.filter(compound_set=compound_set)
    pdb_filepaths = [c.pdb_info.path for c in computed]

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
    # filter_permissions = "target_id__project_id"
    # filter_fields = '__all__'
### End of Session Project
