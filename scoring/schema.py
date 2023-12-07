import graphene
from graphene_django.rest_framework.mutation import SerializerMutation

from scoring.serializers import (
    CmpdChoiceSerializer,
    ScoreChoiceSerializer,
    ViewSceneSerializer,
)

relay = graphene.relay


class ViewScene(SerializerMutation):
    class Meta:
        serializer_class = ViewSceneSerializer
        interfaces = (relay.Node,)


class CmpdChoice(SerializerMutation):
    class Meta:
        serializer_class = CmpdChoiceSerializer
        interfaces = (relay.Node,)


class ScoreChoice(SerializerMutation):
    class Meta:
        serializer_class = ScoreChoiceSerializer
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    view_scene = graphene.List(ViewScene)
    cmpd_choice = graphene.List(CmpdChoice)
    score_choice = graphene.List(ScoreChoice)


schema = graphene.Schema(query=Query)
