from django.http import HttpResponse
from django.shortcuts import render
from viewer.models import Molecule, Protein, Compound
import json
from frag.network.decorate import get_3d_vects_for_mol
from api.utils import draw_mol
from .tasks import get_my_graph
from viewer.models import Molecule, Protein, Compound, Target
from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
)
from rest_framework import viewsets


class ISpyBSafeQuerySet(viewsets.ReadOnlyModelViewSet):

    def get_open_proposals(self):
        """
        Returns the list of proposals anybody can access
        :return:
        """
        return ["PROPOSAL_THREE"]

    def get_proposals_for_user_dummy(self, user):
        DUMMY_DATA = {"qkg34138": "PROPOSAL_ONE", "uzw12877": "PROPOSAL_TWO"}
        return DUMMY_DATA[user]

    def get_proposals_for_user_from_ispyb(self, user):
        import ispyb

        # Check ISPyB for the relavent data
        # Get a connection and data area objects
        # What username and password to use in config.cfg
        with ispyb.open("config.cfg") as conn:
            core = conn.core
            mx_acquisition = conn.mx_acquisition
            return mx_acquisition.get_proposals_for_user()

    def get_proposals_for_user(self):
        user = self.request.user
        return self.get_proposals_for_user_dummy(user)

    def get_queryset(self):
        """
        Optionally restricts the returned purchases to a given propsals
        """
        # The list of proposals this user can have
        proposal_list = self.get_proposals_for_user()
        # Add in the ones everyone has access to
        proposal_list.extend(self.get_open_proposals())
        # Must have a directy foreign key (project_id) for it to work
        return self.queryset.object.filter(project_id__in=proposal_list)


class TargetView(ISpyBSafeQuerySet):
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
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


def get_vects_from_pk(request, pk):
    sdf_info = str(Molecule.objects.get(pk=pk).sdf_info)
    out_data = get_3d_vects_for_mol(sdf_info)
    return HttpResponse(json.dumps(out_data))


def get_graph_from_pk(request, pk):
    smiles = str(Molecule.objects.get(pk=pk).smiles)
    out_data = get_my_graph(smiles)
    return HttpResponse(json.dumps(out_data))
