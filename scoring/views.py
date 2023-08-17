import json
from django.http import HttpResponse
from frag.conf.functions import generate_confs_for_vector
from rest_framework import viewsets

from scoring.models import (
    ViewScene,
    CmpdChoice,
    SiteObservationChoice,
    SiteObservationAnnotation,
    ScoreChoice,
    SiteObservationGroup,
)
from scoring.serializers import (
    ViewSceneSerializer,
    CmpdChoiceSerializer,
    SiteObservationChoiceSerializer,
    SiteObservationAnnotationSerializer,
    ScoreChoiceSerializer,
    SiteObservationGroupSerializer,
)

class ViewSceneView(viewsets.ModelViewSet):
    queryset = ViewScene.objects.filter().order_by('-modified')
    # filter_backends = (filters.DjangoFilterBackend,)
    serializer_class = ViewSceneSerializer
    filter_fields = ("user_id", "uuid")

    def put(self, request, *args, **kwargs):
        return self.partial_update(request, *args, **kwargs)


class SiteObservationChoiceView(viewsets.ModelViewSet):
    queryset = SiteObservationChoice.objects.filter()
    serializer_class = SiteObservationChoiceSerializer
    filter_fields = ("user_id", "mol_id", "mol_id__prot_id__target_id", "choice_type")    


class SiteObservationAnnotationView(viewsets.ModelViewSet):
    queryset = SiteObservationAnnotation.objects.filter()
    serializer_class = SiteObservationAnnotationSerializer
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


class SiteObservationGroupView(viewsets.ModelViewSet):
    queryset = SiteObservationGroup.objects.filter()
    serializer_class = SiteObservationGroupSerializer
    filter_fields = ("group_type", "mol_id", "target_id", "description")    


def gen_conf_from_vect(request):
    input_dict = json.loads(request.body)
    input_smiles = input_dict["INPUT_SMILES"]
    input_mol_block = input_dict["INPUT_MOL_BLOCK"]
    return HttpResponse(
        json.dumps(
            generate_confs_for_vector(input_smiles, input_mol_block)
        )
    )


def get_current_user_id(request):
    if request.user.is_authenticated():
        return HttpResponse(json.dumps(request.user.id))
    else:
        return HttpResponse(json.dumps('null'))
