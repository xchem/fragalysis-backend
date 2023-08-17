from rest_framework import serializers

from scoring.models import (
    ViewScene,
    CmpdChoice,
    SiteObservationChoice,
    ScoreChoice,
    SiteObservationGroup,
    SiteObservationAnnotation,
)


class ViewSceneSerializer(serializers.ModelSerializer):

    class Meta:
        model = ViewScene
        fields = ("id", "uuid", "title", "scene", "created", "modified", "user_id", "snapshot")


class SiteObservationChoiceSerializer(serializers.ModelSerializer):

    class Meta:
        model = SiteObservationChoice
        fields = ("id", "user_id", "mol_id", "choice_type", "score")        


class SiteObservationAnnotationSerializer(serializers.ModelSerializer):

    class Meta:
        model = SiteObservationAnnotation
        fields = ("id", "mol_id", "annotation_type", "annotation_text")        


class CmpdChoiceSerializer(serializers.ModelSerializer):

    class Meta:
        model = CmpdChoice
        fields = ("id", "user_id", "cmpd_id", "choice_type", "score")


class ScoreChoiceSerializer(serializers.ModelSerializer):

    class Meta:
        model = ScoreChoice
        fields = (
            "id",
            "user_id",
            "mol_id",
            "prot_id",
            "choice_type",
            "score",
            "is_done",
        )


class SiteObservationGroupSerializer(serializers.ModelSerializer):

    class Meta:
        model = SiteObservationGroup
        fields = (
            "id",
            "group_type",
            "mol_id",
            "target_id",
            "x_com",
            "y_com",
            "z_com",
            "description",
        )        


