import json, os

from django.db import connections
from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.conf import settings

from api.security import ISpyBSafeQuerySet
from api.utils import get_params, get_highlighted_diffs
from viewer.models import Molecule, Protein, Compound, Target, SessionProject, Snapshot, ComputedCompound
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
    choice = 'none yet'
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

            name = myfile.name
            path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
            tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

            # isfile = os.path.isfile(tmp_file)
            d, v = validate(tmp_file, target=target)
            print(d)
            print(v)
            pd.set_option('display.max_colwidth', -1)
            if not v:
                table = pd.DataFrame.from_dict(d)
                html_table = table.to_html()
                html_table += '''<p> Your data was <b>not</b> validated. The table above shows errors</p>'''
                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html_table})
                # return ValidationError('We could not validate this file')
            if choice=='upload':
                cset = process_compound_set(target=target, filename=tmp_file)
                os.remove(tmp_file)

                computed = ComputedCompound.objects.filter(compound_set=cset).values()
                table = pd.DataFrame(computed)
                html_table = table.to_html()
                html_table += '''<p> Your data was validated. The table above shows the compounds in the set</p>'''

                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html_table})
            if choice=='validate' and v:
                html = '<p> Your data was validated. You can upload it by checking the upload radio button</p>'
                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html})
    else:
        # GET, generate blank form
        form = CSetForm()
        # return render(request, 'viewer/upload-cset.html', {
        #     'uploaded_file_url': uploaded_file_url
        # })
    return render(request, 'viewer/upload-cset.html', {'form': form, 'table': choice})


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
    filter_permissions = "target_id__project_id"
    filter_fields = '__all__'
### End of Session Project

# from rest_framework.parsers import FileUploadParser
# from rest_framework.response import Response
# from rest_framework.views import APIView
# from rest_framework import status
#
# from .serializers import FileSerializer
#
#
# class FileUploadView(APIView):
#     parser_class = (FileUploadParser,)
#
#     def post(self, request, *args, **kwargs):
#
#       file_serializer = FileSerializer(data=request.data)
#
#       if file_serializer.is_valid():
#           file_serializer.save()
#           return Response(file_serializer.data, status=status.HTTP_201_CREATED)
#       else:
#           return Response(file_serializer.errors, status=status.HTTP_400_BAD_REQUEST)