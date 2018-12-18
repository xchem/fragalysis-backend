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
    paginate_by = None
    queryset = ViewScene.objects.filter().order_by('-modified')
    serializer_class = ViewSceneSerializer
    filter_fields = ("user_id", "uuid")

    def put(self, request, *args, **kwargs):
        return self.partial_update(request, *args, **kwargs)


class ProtChoiceView(viewsets.ModelViewSet):
    paginate_by = None
    queryset = ProtChoice.objects.filter()
    serializer_class = ProtChoiceSerializer
    filter_fields = ("user_id", "prot_id", "prot_id__target_id", "choice_type")


class MolChoiceView(viewsets.ModelViewSet):
    paginate_by = None
    queryset = MolChoice.objects.filter()
    serializer_class = MolChoiceSerializer
    filter_fields = ("user_id", "mol_id", "mol_id__prot_id__target_id", "choice_type")


class MolAnnotationView(viewsets.ModelViewSet):
    paginate_by = None
    queryset = MolAnnotation.objects.filter()
    serializer_class = MolAnnotationSerializer
    filter_fields = ("mol_id", "annotation_type")


class CmpdChoiceView(viewsets.ModelViewSet):
    paginate_by = None
    queryset = CmpdChoice.objects.filter()
    serializer_class = CmpdChoiceSerializer
    filter_fields = ("user_id", "cmpd_id", "choice_type")


class ScoreChoiceView(viewsets.ModelViewSet):
    paginate_by = None
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
    paginate_by = None
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


def get_current_user_id(request):
    if request.user.is_authenticated():
        return HttpResponse(json.dumps(request.user.id))
    else:
        return HttpResponse(json.dumps('null'))
