from rest_framework import serializers
from hypothesis.models import (
    Vector3D,
    Vector,
    Interaction,
    InteractionPoint,
    ProteinResidue,
    TargetResidue,
)


class Vector3DSerializer(serializers.ModelSerializer):

    class Meta:
        model = Vector3D
        fields = (
            "id",
            "mol_id",
            "vector_id",
            "number",
            "start_x",
            "start_y",
            "start_z",
            "end_x",
            "end_y",
            "end_z",
        )


class VectorSerializer(serializers.ModelSerializer):

    class Meta:
        model = Vector
        fields = ("id", "cmpd_id", "smiles", "type")


class InteractionSerializer(serializers.ModelSerializer):

    class Meta:
        model = Interaction
        fields = (
            "id",
            "interaction_version",
            "interaction_point",
            "interaction_type",
            "distance",
            "score",
            "prot_smarts",
            "mol_smarts",
        )


class InteractionPointSerializer(serializers.ModelSerializer):

    class Meta:
        model = InteractionPoint
        fields = (
            "id",
            "prot_res_id",
            "mol_id",
            "protein_atom_name",
            "molecule_atom_name",
        )


class ProteinResidueSerialzier(serializers.ModelSerializer):

    class Meta:
        model = ProteinResidue
        fields = ("id", "prot_id", "targ_res_id")


class TargetResidueSerialzier(serializers.ModelSerializer):

    class Meta:
        model = TargetResidue
        fields = ("id", "target_id", "res_name", "res_num", "chain_id")
