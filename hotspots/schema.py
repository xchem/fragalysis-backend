import graphene
from graphene_django.rest_framework.mutation import SerializerMutation
from hotspots.serializers import HotspotMapSerializer

relay = graphene.relay


class HotspotMap(SerializerMutation):

    class Meta:
        serializer_class = HotspotMapSerializer
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    hotspotmap = graphene.List(HotspotMap)


schema = graphene.Schema(query=Query)
