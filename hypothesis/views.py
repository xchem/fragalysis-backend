from rest_framework import viewsets

from hypothesis.models import (
    Vector,
    Vector3D,
    Interaction,
    InteractionPoint,
    TargetResidue,
)
from hypothesis.serializers import (
    VectorSerializer,
    Vector3DSerializer,
    InteractionSerializer,
    InteractionPointSerializer,
    TargetResidueSerialzier,
)


class VectorView(viewsets.ReadOnlyModelViewSet):
    queryset = Vector.objects.filter()
    serializer_class = VectorSerializer
    filterset_fields = ("cmpd_id", "smiles", "type")


class Vector3DView(viewsets.ReadOnlyModelViewSet):
    queryset = Vector3D.objects.filter()
    serializer_class = Vector3DSerializer
    filterset_fields = ("mol_id", "vector_id", "number")


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
