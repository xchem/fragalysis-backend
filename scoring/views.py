import json

from django.http import HttpResponse
from frag.conf.functions import generate_confs_for_vector
from rest_framework import viewsets

from scoring.models import (
    CmpdChoice,
    ScoreChoice,
    SiteObservationAnnotation,
    SiteObservationChoice,
    SiteObservationGroup,
    ViewScene,
)
from scoring.serializers import (
    CmpdChoiceSerializer,
    ScoreChoiceSerializer,
    SiteObservationAnnotationSerializer,
    SiteObservationChoiceSerializer,
    SiteObservationGroupSerializer,
    ViewSceneSerializer,
)


class ViewSceneView(viewsets.ModelViewSet):
    queryset = ViewScene.objects.all().order_by("-modified")
    # filter_backends = (filters.DjangoFilterBackend,)
    serializer_class = ViewSceneSerializer
    filterset_fields = ("user_id", "uuid")

    def put(self, request, *args, **kwargs):
        return self.partial_update(request, *args, **kwargs)


class SiteObservationChoiceView(viewsets.ModelViewSet):
    queryset = SiteObservationChoice.objects.all()
    serializer_class = SiteObservationChoiceSerializer
    filterset_fields = (
        "user",
        "site_observation",
        "site_observation__experiment__experiment_upload__target",
        "choice_type",
    )


class SiteObservationAnnotationView(viewsets.ModelViewSet):
    queryset = SiteObservationAnnotation.objects.all()
    serializer_class = SiteObservationAnnotationSerializer
    filterset_fields = ("site_observation", "annotation_type")


class CmpdChoiceView(viewsets.ModelViewSet):
    queryset = CmpdChoice.objects.all()
    serializer_class = CmpdChoiceSerializer
    filterset_fields = ("user_id", "cmpd_id", "choice_type")


class ScoreChoiceView(viewsets.ModelViewSet):
    queryset = ScoreChoice.filter_manager.filter_qs()
    serializer_class = ScoreChoiceSerializer
    filterset_fields = (
        "user",
        "site_observation",
        "is_done",
        "site_observation__experiment__experiment_upload__target",
        "choice_type",
    )


class SiteObservationGroupView(viewsets.ModelViewSet):
    queryset = SiteObservationGroup.objects.all()
    serializer_class = SiteObservationGroupSerializer
    filterset_fields = ("group_type", "site_observation", "target", "description")


def gen_conf_from_vect(request):
    input_dict = json.loads(request.body)
    input_smiles = input_dict["INPUT_SMILES"]
    input_mol_block = input_dict["INPUT_MOL_BLOCK"]
    return HttpResponse(
        json.dumps(generate_confs_for_vector(input_smiles, input_mol_block))
    )


def get_current_user_id(request):
    if request.user.is_authenticated():
        return HttpResponse(json.dumps(request.user.id))
    else:
        return HttpResponse(json.dumps("null"))
