import json

from django.db import connections
from django.http import HttpResponse
from django.shortcuts import render

from api.security import ISpyBSafeQuerySet
from api.utils import get_params
from viewer.models import Molecule, Protein, Compound, Target
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
)


class VectorsView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Molecule.objects.filter()
    serializer_class = VectorsSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class GraphView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Molecule.objects.filter()
    serializer_class = GraphSerializer
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class MolImageView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Molecule.objects.filter()
    serializer_class = MolImageSerialzier
    filter_permissions = "prot_id__target_id__project_id"
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class CompoundImageView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Compound.objects.filter()
    serializer_class = CmpdImageSerialzier
    filter_permissions = "project_id"
    filter_fields = ("smiles",)


class ProteinMapInfoView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Protein.objects.filter()
    serializer_class = ProtMapInfoSerialzer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


class ProteinPDBInfoView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Protein.objects.filter()
    serializer_class = ProtPDBInfoSerialzer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


class ProteinPDBBoundInfoView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Protein.objects.filter()
    serializer_class = ProtPDBBoundInfoSerialzer
    filter_permissions = "target_id__project_id"
    filter_fields = ("code", "target_id", "prot_type")


class TargetView(ISpyBSafeQuerySet):
    paginate_by = None
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    filter_permissions = "project_id"
    filter_fields = ("title",)


class MoleculeView(ISpyBSafeQuerySet):
    paginate_by = None
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
    paginate_by = None
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_permissions = "project_id"
    filter_fields = ("smiles",)


class ProteinView(ISpyBSafeQuerySet):
    paginate_by = None
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


def img_from_smiles(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        if smiles:
            return get_params(smiles, request)
        else:
            return HttpResponse("Please insert SMILES")
    else:
        return HttpResponse("Please insert SMILES")


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
