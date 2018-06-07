from viewer.models import Molecule, Protein, Compound, Target
from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
)
from rest_framework import permissions
from rest_framework import viewsets


class TargetView(viewsets.ReadOnlyModelViewSet):
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    # permission_classes =  [permissions.DjangoObjectPermissions,]
    filter_fields = ("title",)


class MoleculeView(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    filter_fields = ("prot_id", "cmpd_id", "smiles", "prot_id__target_id", "mol_groups")


class CompoundView(viewsets.ReadOnlyModelViewSet):
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_fields = ("smiles",)


class ProteinView(viewsets.ReadOnlyModelViewSet):
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    filter_fields = ("code", "target_id")


from django.http import HttpResponse
from django.shortcuts import render
from viewer.models import Molecule, Protein, Compound
from scoring.models import MolChoice, ViewScene
from uuid import uuid4
import json
from frag.network.decorate import get_3d_vects_for_mol
from api.utils import draw_mol, get_token
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from .tasks import get_my_graph


def react(request):
    """
    :param request:
    :return:
    """
    return render(request, "viewer/react_temp.html")


def get_params(smiles, request):
    height = None
    if "height" in request.GET:
        height = int(request.GET["height"])
    width = None
    if "width" in request.GET:
        width = int(request.GET["width"])
    return HttpResponse(draw_mol(smiles, width=width, height=height))


def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"].rstrip(".svg")
        return get_params(smiles, request)
    else:
        return HttpResponse("Please insert SMILES")


def img_from_pk(request):
    """
    DEPRECATED DESTROY!!!
    :param request:
    :return:
    """
    if "pk" in request.GET:
        smiles = Molecule.objects.get(pk=request.GET["pk"]).smiles
        return get_params(smiles, request)
    else:
        return HttpResponse("Please insert PK")


def img_from_mol_pk(request, pk):
    smiles = Molecule.objects.get(pk=pk).smiles
    return get_params(smiles, request)


def img_from_cmpd_pk(request, pk):
    smiles = Compound.objects.get(pk=pk).smiles
    return get_params(smiles, request)


def img_from_smiles(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        return get_params(smiles, request)
    else:
        return HttpResponse("Please insert SMILES")


def mol_from_pk(request, pk):
    sdf_info = Molecule.objects.get(pk=pk).sdf_info
    return HttpResponse(sdf_info)


def map_from_pk(request, pk):
    prot = Protein.objects.get(pk=pk)
    return HttpResponse(prot.map_info)


def prot_from_pk(request, pk):
    pdb_info = open(Protein.objects.get(pk=pk).pdb_info.path).read()
    return HttpResponse(pdb_info)


def post_view(request):
    """
    Post the view for a given scene
    :param request:
    :return:
    """
    new_view = ViewScene()
    new_view.title = request.POST["title"]
    new_view.scene = request.POST["scene"]
    new_view.uuid = str(uuid4())
    new_view.save()
    url = request.build_absolute_uri("?scene_id=" + str(new_view.pk)).replace(
        "/post_view/", "/display/"
    )
    return HttpResponse(url)


def get_view(request, pk):
    """
    Now Get the view for a given UUID
    :param request:
    :return:
    """
    this_view = ViewScene.objects.get(pk=pk)
    return HttpResponse(
        json.dumps({"title": this_view.title, "scene": this_view.scene})
    )


def get_vects_from_pk(request, pk):
    sdf_info = str(Molecule.objects.get(pk=pk).sdf_info)
    out_data = get_3d_vects_for_mol(sdf_info)
    return HttpResponse(json.dumps(out_data))


def get_graph_from_pk(request, pk):
    smiles = str(Molecule.objects.get(pk=pk).smiles)
    out_data = get_my_graph(smiles)
    return HttpResponse(json.dumps(out_data))
