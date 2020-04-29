import json, os
import zipfile

from django.db import connections
from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.conf import settings
from django.contrib.sites.shortcuts import get_current_site

from api.security import ISpyBSafeQuerySet
from api.utils import get_params, get_highlighted_diffs

from viewer.models import Molecule, Protein, Compound, Target, SessionProject, Snapshot, ComputedCompound, CompoundSet
from viewer import filters
from sdf_check import validate
from forms import CSetForm
from compound_sets import process_compound_set
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


# needs a target to be specified
def upload_cset(request):
    """
    :param request:
    :return:
    """
    zfile = None
    zf = None
    if request.method == 'POST':
        # POST, generate form with data from the request
        print('data provided... processing')
        form = CSetForm(request.POST, request.FILES)
        # check if it's valid:
        if form.is_valid():
            myfile = request.FILES['sdf_file']
            print(myfile)
            target = request.POST['target_name']
            choice = request.POST['submit_choice']
            if 'pdb_zip' in request.FILES.keys():
                pdb_file = request.FILES['pdb_zip']
            else:
                pdb_file = None

            # if request.FILES['pdb_zip']!='':
                # check it's actually a zip file

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


            name = myfile.name
            path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
            tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

            # isfile = os.path.isfile(tmp_file)
            d, v = validate(tmp_file, target=target, zfile=zfile)
            print(d)
            print(v)
            pd.set_option('display.max_colwidth', -1)
            if not v:
                table = pd.DataFrame.from_dict(d)
                html_table = table.to_html()
                html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''
                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html_table})
                # return ValidationError('We could not validate this file')
            if str(choice)=='1':
                cset = process_compound_set(target=target, filename=tmp_file, zfile=zfile)
                if zf:
                    zf.close()
                # computed = ComputedCompound.objects.filter(compound_set=cset).values()
                submitter = cset.submitter
                name = submitter.unique_name

                download_url = '<a href="/viewer/compound_set/%s">Download Compound Set</a>' %name

                # table = pd.DataFrame(computed)
                # html_table = table.to_html()
                html_table = '''<p> Your data was validated and added to the fragalysis database. The link above will allow you to download the submitted file</p>'''

                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html_table, 'download_url': download_url})
            if str(choice)=='0' and v:
                html = '<p> Your data was validated. You can upload it by checking the upload radio button</p>'
                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html, 'download_url':''})
    else:
        # GET, generate blank form
        form = CSetForm()
        # return render(request, 'viewer/upload-cset.html', {
        #     'uploaded_file_url': uploaded_file_url
        # })
    return render(request, 'viewer/upload-cset.html', {'form': form, 'table': '', 'download_url':''})


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
