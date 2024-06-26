from rest_framework import viewsets

from api.security import ISPyBSafeQuerySet
from hypothesis.models import Interaction, InteractionPoint, TargetResidue
from hypothesis.serializers import (
    InteractionPointSerializer,
    InteractionSerializer,
    TargetResidueSerialzier,
)


class InteractionView(viewsets.ReadOnlyModelViewSet):
    queryset = Interaction.objects.filter()
    serializer_class = InteractionSerializer
    filterset_fields = (
        "interaction_point",
        "interaction_version",
        "interaction_type",
        "interaction_point__targ_res__target_id",
        "interaction_point__site_observation",
        "distance",
        "score",
        "prot_smarts",
        "mol_smarts",
    )


class InteractionPointView(ISPyBSafeQuerySet):
    queryset = InteractionPoint.objects.all()
    serializer_class = InteractionPointSerializer
    filterset_fields = ("site_observation", "protein_atom_name", "molecule_atom_name")
    filter_permissions = "targ_res__target_id__project_id"


class TargetResidueView(ISPyBSafeQuerySet):
    queryset = TargetResidue.objects.all()
    serializer_class = TargetResidueSerialzier
    filterset_fields = ("target_id", "res_name", "res_num", "chain_id")
    filter_permissions = "target_id__project_id"
