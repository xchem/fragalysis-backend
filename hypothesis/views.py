
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
        "interaction_point__prot_res_id__targ_res_id__target_id",
        "interaction_point__mol_id",
        "interaction_point__prot_res_id__prot_id",
        "distance",
        "score",
        "prot_smarts",
        "mol_smarts",
    )


class InteractionPointView(viewsets.ReadOnlyModelViewSet):
    queryset = InteractionPoint.objects.filter()
    serializer_class = InteractionPointSerializer
    filterset_fields = ("prot_res_id", "mol_id", "protein_atom_name", "molecule_atom_name")


class TargetResidueView(viewsets.ReadOnlyModelViewSet):
    queryset = TargetResidue.objects.filter()
    serializer_class = TargetResidueSerialzier
    filterset_fields = ("target_id", "res_name", "res_num", "chain_id")
