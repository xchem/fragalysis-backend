from frag.conf.functions import generate_confs_for_vector
from django.http import HttpResponse
import json
from scoring.models import (
    ViewScene,
    ProtChoice,
    CmpdChoice,
    MolChoice,
    ScoreChoice,
    MolGroup,
)
from scoring.serializers import (
    ViewSceneSerializer,
    ProtChoiceSerializer,
    CmpdChoiceSerializer,
    MolChoiceSerializer,
    ScoreChoiceSerializer,
    MolGroupSerializer,
)
from rest_framework import viewsets


class ViewSceneView(viewsets.ModelViewSet):
    queryset = ViewScene.objects.filter()
    serializer_class = ViewSceneSerializer
    filter_fields = ("user_id", "uuid")


class ProtChoiceView(viewsets.ModelViewSet):
    queryset = ProtChoice.objects.filter()
    serializer_class = ProtChoiceSerializer
    filter_fields = ("user_id", "prot_id", "prot_id__target_id", "choice_type")


class MolChoiceView(viewsets.ModelViewSet):
    queryset = MolChoice.objects.filter()
    serializer_class = MolChoiceSerializer
    filter_fields = ("user_id", "mol_id", "mol_id__prot_id__target_id", "choice_type")


class CmpdChoiceView(viewsets.ModelViewSet):
    queryset = CmpdChoice.objects.filter()
    serializer_class = CmpdChoiceSerializer
    filter_fields = ("user_id", "cmpd_id", "choice_type")


class ScoreChoiceView(viewsets.ModelViewSet):
    queryset = ScoreChoice.objects.filter()
    serializer_class = ScoreChoiceSerializer
    filter_fields = (
        "user_id",
        "mol_id",
        "prot_id",
        "is_done",
        "mol_id__prot_id__target_id",
        "prot_id__target_id",
        "choice_type",
    )


class MolGroupView(viewsets.ModelViewSet):
    queryset = MolGroup.objects.filter()
    serializer_class = MolGroupSerializer
    filter_fields = ("group_type", "mol_id", "target_id")


def gen_conf_from_vect(request):
    input_vector = request.POST["INPUT_VECTOR"]
    input_smiles = request.POST.getlist("INPUT_SMILES")
    input_mol_block = request.POST["INPUT_MOL_BLOCK"]
    return HttpResponse(json.dumps(generate_confs_for_vector(input_vector, input_smiles, input_mol_block)))
