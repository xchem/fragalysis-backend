import json

from django.http import HttpResponse
from frag.conf.functions import generate_confs_for_vector
from rest_framework import mixins

from api.security import ISPyBSafeQuerySet
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
from viewer.permissions import IsObjectProposalMember


class ViewSceneView(
    mixins.UpdateModelMixin,
    mixins.CreateModelMixin,
    mixins.DestroyModelMixin,
    ISPyBSafeQuerySet,
):
    queryset = ViewScene.objects.all().order_by("-modified")
    # filter_backends = (filters.DjangoFilterBackend,)
    serializer_class = ViewSceneSerializer
    filterset_fields = ("user_id", "uuid")
    filter_permissions = "snapshot__session_project__target__project_id"
    permission_classes = [IsObjectProposalMember]

    def put(self, request, *args, **kwargs):
        return self.partial_update(request, *args, **kwargs)


class SiteObservationChoiceView(
    mixins.UpdateModelMixin,
    mixins.CreateModelMixin,
    mixins.DestroyModelMixin,
    ISPyBSafeQuerySet,
):
    queryset = SiteObservationChoice.objects.all()
    serializer_class = SiteObservationChoiceSerializer
    filterset_fields = (
        "user",
        "site_observation",
        "site_observation__experiment__experiment_upload__target",
        "choice_type",
    )
    filter_permissions = "site_observation__experiment__experiment_upload__project"
    permission_classes = [IsObjectProposalMember]


class SiteObservationAnnotationView(
    mixins.UpdateModelMixin,
    mixins.CreateModelMixin,
    mixins.DestroyModelMixin,
    ISPyBSafeQuerySet,
):
    queryset = SiteObservationAnnotation.objects.all()
    serializer_class = SiteObservationAnnotationSerializer
    filterset_fields = ("site_observation", "annotation_type")
    filter_permissions = "site_observation__experiment__experiment_upload__project"
    permission_classes = [IsObjectProposalMember]


class CmpdChoiceView(
    mixins.UpdateModelMixin,
    mixins.CreateModelMixin,
    mixins.DestroyModelMixin,
    ISPyBSafeQuerySet,
):
    queryset = CmpdChoice.objects.all()
    serializer_class = CmpdChoiceSerializer
    filterset_fields = ("user_id", "cmpd_id", "choice_type")
    filter_permissions = "cmpd_id__project_id"
    permission_classes = [IsObjectProposalMember]


class ScoreChoiceView(
    mixins.UpdateModelMixin,
    mixins.CreateModelMixin,
    mixins.DestroyModelMixin,
    ISPyBSafeQuerySet,
):
    queryset = ScoreChoice.filter_manager.filter_qs()
    serializer_class = ScoreChoiceSerializer
    filterset_fields = (
        "user",
        "site_observation",
        "is_done",
        "site_observation__experiment__experiment_upload__target",
        "choice_type",
    )
    filter_permissions = "site_observation__experiment__experiment_upload__project"
    permission_classes = [IsObjectProposalMember]


class SiteObservationGroupView(
    mixins.UpdateModelMixin,
    mixins.CreateModelMixin,
    mixins.DestroyModelMixin,
    ISPyBSafeQuerySet,
):
    queryset = SiteObservationGroup.objects.all()
    serializer_class = SiteObservationGroupSerializer
    filterset_fields = ("group_type", "site_observation", "target", "description")
    filter_permissions = "target__project_id"
    permission_classes = [IsObjectProposalMember]


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
