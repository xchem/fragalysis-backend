from rest_framework import viewsets
from rest_framework.views import APIView
from rest_framework.response import Response

from hypothesis.models import (
    Interaction,
    InteractionPoint,
    TargetResidue,
)
from hypothesis.serializers import (
    InteractionSerializer,
    InteractionPointSerializer,
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


class InteractionPointView(viewsets.ReadOnlyModelViewSet):
    queryset = InteractionPoint.objects.all()
    serializer_class = InteractionPointSerializer
    filterset_fields = ("site_observation", "protein_atom_name", "molecule_atom_name")


class TargetResidueView(viewsets.ReadOnlyModelViewSet):
    queryset = TargetResidue.objects.all()
    serializer_class = TargetResidueSerialzier
    filterset_fields = ("target_id", "res_name", "res_num", "chain_id")
