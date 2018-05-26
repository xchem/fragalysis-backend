from graphene_django import DjangoObjectType
import graphene
from graphene import relay
from graphene_django.rest_framework.mutation import SerializerMutation
from pandda.serializers import PanddaEventSerializer, PanddaSiteSerializer


class PanddaEvent(SerializerMutation):

    class Meta:
        serializer_class = PanddaEventSerializer
        interfaces = (relay.Node,)


class PanddaSite(SerializerMutation):

    class Meta:
        serializer_class = PanddaSiteSerializer
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    pandda_event = graphene.List(PanddaEvent)
    pandda_site = graphene.List(PanddaSite)


schema = graphene.Schema(query=Query)
