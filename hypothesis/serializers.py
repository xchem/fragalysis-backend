from rest_framework import serializers

from hypothesis.models import (
    Vector3D,
    Vector,
    Interaction,
    InteractionPoint,
    TargetResidue,
)
from viewer.models import SiteObservation


# class Vector3DSerializer(serializers.ModelSerializer):

#     class Meta:
#         model = Vector3D
#         fields = (
#             "id",
#             "site_observation",
#             "vector_id",
#             "number",
#             "start_x",
#             "start_y",
#             "start_z",
#             "end_x",
#             "end_y",
#             "end_z",
#         )

class Vector3DSerializer(serializers.Serializer):

    start_x = serializers.FloatField()
    start_y = serializers.FloatField()
    start_z = serializers.FloatField()
    end_x = serializers.FloatField()
    end_y = serializers.FloatField()
    end_z = serializers.FloatField()
    number = serializers.IntegerField()
    vector_type = serializers.CharField()
    smiles = serializers.CharField()
    site_observation = serializers.IntegerField()
    cmpd_id = serializers.IntegerField()
    


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


class TargetResidueSerialzier(serializers.ModelSerializer):

    class Meta:
        model = TargetResidue
        fields = ("id", "target_id", "res_name", "res_num", "chain_id")
