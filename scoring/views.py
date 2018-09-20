import json

from django.http import HttpResponse
from frag.conf.functions import generate_confs_for_vector
from rest_framework import viewsets

from scoring.models import (
    ViewScene,
    ProtChoice,
    CmpdChoice,
    MolChoice,
    MolAnnotation,
    ScoreChoice,
    MolGroup,
)
from scoring.serializers import (
    ViewSceneSerializer,
    ProtChoiceSerializer,
    CmpdChoiceSerializer,
    MolChoiceSerializer,
    MolAnnotationSerializer,
    ScoreChoiceSerializer,
    MolGroupSerializer,
)


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


class MolAnnotationView(viewsets.ModelViewSet):
    queryset = MolAnnotation.objects.filter()
    serializer_class = MolAnnotationSerializer
    filter_fields = ("mol_id", "annotation_type")


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
    filter_fields = ("group_type", "mol_id", "target_id", "description")


def gen_conf_from_vect(request):
    input_dict = json.loads(request.body)
    input_vector = sorted(input_dict["INPUT_VECTOR"].split("."), reverse=True)[
        0
    ].replace("Xe", "H")
    input_smiles = input_dict["INPUT_SMILES"]
    input_mol_block = input_dict["INPUT_MOL_BLOCK"]
    return HttpResponse(
        json.dumps(
            generate_confs_for_vector(input_vector, input_smiles, input_mol_block)
        )
    )
