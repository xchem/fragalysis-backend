import graphene
from graphene import relay
from graphene_django.rest_framework.mutation import SerializerMutation

from scoring.serializers import (
    ViewSceneSerializer,
    ProtChoiceSerializer,
    CmpdChoiceSerializer,
    MolChoiceSerializer,
    ScoreChoiceSerializer,
    MolGroupSerializer,
)


class ViewScene(SerializerMutation):

    class Meta:
        serializer_class = ViewSceneSerializer
        interfaces = (relay.Node,)


class ProtChoice(SerializerMutation):

    class Meta:
        serializer_class = ProtChoiceSerializer
        interfaces = (relay.Node,)


class MolChoice(SerializerMutation):

    class Meta:
        serializer_class = MolChoiceSerializer
        interfaces = (relay.Node,)


class CmpdChoice(SerializerMutation):

    class Meta:
        serializer_class = CmpdChoiceSerializer
        interfaces = (relay.Node,)


class ScoreChoice(SerializerMutation):

    class Meta:
        serializer_class = ScoreChoiceSerializer
        interfaces = (relay.Node,)


class MolGroup(SerializerMutation):

    class Meta:
        serializer_class = MolGroupSerializer
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    view_scene = graphene.List(ViewScene)
    prot_choice = graphene.List(ProtChoice)
    mol_choice = graphene.List(MolChoice)
    cmpd_choice = graphene.List(CmpdChoice)
    score_choice = graphene.List(ScoreChoice)
    mol_group = graphene.List(MolGroup)


schema = graphene.Schema(query=Query)
