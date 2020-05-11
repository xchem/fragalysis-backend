import json, os
import zipfile
from cStringIO import StringIO

from django.db import connections
from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets, views
from rest_framework.response import Response
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.conf import settings

from api.security import ISpyBSafeQuerySet
from api.utils import get_params, get_highlighted_diffs

from viewer.models import (
    Molecule,
    Protein,
    Compound,
    Target,
    SessionProject,
    Snapshot,
    ComputedCompound,
    CompoundSet,
    CSetKeys,
    NumericalScoreValues,
    ScoreDescription
)
from viewer import filters
from sdf_check import validate
from forms import CSetForm, UploadKeyForm
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
    SnapshotWriteSerializer,
    TargetCompoundSetsSerializer,
    CompoundSetSerializer,
    CompoundMoleculeSerializer,
    NumericalScoreSerializer,
    ScoreDescriptionAllSerializer
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


# needs a target to be specified
def upload_cset(request):
    """
    :param request:
    :return:
    """
    zfile = None
    zf = None
    cset = None
    if request.method == 'POST':
        form = CSetForm(request.POST, request.FILES)
        # POST, generate form with data from the request
        key = request.POST['upload_key']
        all_keys = CSetKeys.objects.all()
        if key not in [str(key.uuid) for key in all_keys]:
            html = "<br><p>You either didn't provide an upload key, or it wasn't valid. Please try again (email rachael.skyner@diamond.ac.uk to obtain an upload key)</p>"
            return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html, 'download_url':''})
        print('data provided... processing')
#         try:
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
                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html_table, 'download_url':''})
                # return ValidationError('We could not validate this file')
            if str(choice)=='1':
                cset = process_compound_set(target=target, filename=tmp_file, zfile=zfile)
                if zf:
                    zf.close()
                # computed = ComputedCompound.objects.filter(compound_set=cset).values()
                submitter = cset.submitter
                name = submitter.unique_name

                cset_download_url = '<a href="/viewer/compound_set/%s">Download Compound Set</a>' %name
                pset_download_url = '<a href="/viewer/protein_set/%s">Download Protein Set</a>' % name

                download_url = cset_download_url + '<br>' + pset_download_url

                # table = pd.DataFrame(computed)
                # html_table = table.to_html()
                html_table = '''<p> Your data was validated and added to the fragalysis database. The link above will allow you to download the submitted file</p>'''

                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html_table, 'download_url': download_url})
            if str(choice)=='0' and v:
                html = '<p> Your data was validated. You can upload it by checking the upload radio button</p>'
                return render(request, 'viewer/upload-cset.html', {'form': form, 'table': html, 'download_url':''})
#         except:
#             if cset:
#                 cset.delete()
#                 computed = ComputedCompound.objects.filter(compound_set=cset)
#                 for c in computed:
#                     c.delete()
#             return HttpResponse(status=500)
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
### End of Session Project


class TargetCompoundSetsView(views.APIView):
    def get(self, request, targetID):
        target = Target.objects.get(id=targetID)
        compound_sets = CompoundSet.objects.filter(target=target)
        a_compound_set = compound_sets[0]
        compound_mols = ComputedCompound.objects.filter(compound_set=a_compound_set)
        numerical_scores = NumericalScoreValues.objects.filter(compound=compound_mols[0])

        results = TargetCompoundSetsSerializer(numerical_scores, many=True).data
        return Response(results)

class CompoundSetView(viewsets.ReadOnlyModelViewSet):
    queryset = CompoundSet.objects.filter()
    serializer_class = CompoundSetSerializer
    filter_permissions = "project_id"
    filter_fields = ('target',)


class CompoundMoleculesView(viewsets.ReadOnlyModelViewSet):
    queryset = ComputedCompound.objects.filter()
    serializer_class = CompoundMoleculeSerializer
    filter_permissions = "project_id"
    filter_fields = ('compound_set',)

class NumericalScoresView(viewsets.ReadOnlyModelViewSet):
    queryset = NumericalScoreValues.objects.filter()
    serializer_class = NumericalScoreSerializer
    filter_permissions = "project_id"
    filter_fields = ('compound',)


class CompoundScoresView(viewsets.ReadOnlyModelViewSet):
    queryset = ScoreDescription.objects.filter()
    serializer_class = ScoreDescriptionAllSerializer
    filter_permissions = "project_id"
    filter_fields = ('compound_set',)