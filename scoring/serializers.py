from rest_framework import serializers

from scoring.models import (
    CmpdChoice,
    ScoreChoice,
    SiteObservationAnnotation,
    SiteObservationChoice,
    SiteObservationGroup,
    ViewScene,
)


class ViewSceneSerializer(serializers.ModelSerializer):
    class Meta:
        model = ViewScene
        fields = (
            "id",
            "uuid",
            "title",
            "scene",
            "created",
            "modified",
            "user_id",
            "snapshot",
        )


class SiteObservationChoiceSerializer(serializers.ModelSerializer):
    class Meta:
        model = SiteObservationChoice
        fields = ("id", "user", "site_observation", "choice_type", "score")


class SiteObservationAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = SiteObservationAnnotation
        fields = ("id", "site_observation", "annotation_type", "annotation_text")


class CmpdChoiceSerializer(serializers.ModelSerializer):
    class Meta:
        model = CmpdChoice
        fields = ("id", "user_id", "cmpd_id", "choice_type", "score")


class ScoreChoiceSerializer(serializers.ModelSerializer):
    class Meta:
        model = ScoreChoice
        fields = (
            "id",
            "user",
            "site_observation",
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
            "site_observation",
            "target",
            "x_com",
            "y_com",
            "z_com",
            "description",
        )
