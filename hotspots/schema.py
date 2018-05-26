from graphene_django import DjangoObjectType
import graphene
from graphene_django.rest_framework.mutation import SerializerMutation
from hotspots.serializers import HotspotMapSerializer


class HotspotMap(SerializerMutation):

    class Meta:
        serializer_class = HotspotMapSerializer


class Query(graphene.ObjectType):
    hotspotmap = graphene.List(HotspotMap)


schema = graphene.Schema(query=Query)
