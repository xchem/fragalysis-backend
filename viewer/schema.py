import graphene
from graphene_django.rest_framework.mutation import SerializerMutation

from viewer.serializers import CompoundSerializer, TargetSerializer

relay = graphene.relay


class Compound(SerializerMutation):
    class Meta:
        serializer_class = CompoundSerializer
        interfaces = (relay.Node,)


class Target(SerializerMutation):
    class Meta:
        serializer_class = TargetSerializer
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    compounds = graphene.List(Compound)
    targets = graphene.List(Target)


schema = graphene.Schema(query=Query)
