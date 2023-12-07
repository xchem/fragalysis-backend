import graphene
from graphene_django.rest_framework.mutation import SerializerMutation

from hypothesis.serializers import (
    InteractionPointSerializer,
    InteractionSerializer,
    TargetResidueSerialzier,
    Vector3DSerializer,
    VectorSerializer,
)

relay = graphene.relay


class Vector3D(SerializerMutation):
    class Meta:
        serializer_class = Vector3DSerializer
        interfaces = (relay.Node,)


class Vector(SerializerMutation):
    class Meta:
        serializer_class = VectorSerializer
        interfaces = (relay.Node,)


class InteractionPoint(SerializerMutation):
    class Meta:
        serializer_class = InteractionPointSerializer
        interfaces = (relay.Node,)


class Interaction(SerializerMutation):
    class Meta:
        serializer_class = InteractionSerializer
        interfaces = (relay.Node,)


class TargetResidue(SerializerMutation):
    class Meta:
        serializer_class = TargetResidueSerialzier
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    vector3d = graphene.List(Vector3D)
    vector = graphene.List(Vector)
    interaction_point = graphene.List(InteractionPoint)
    interaction = graphene.List(Interaction)
    target_residue = graphene.List(TargetResidue)


schema = graphene.Schema(query=Query)
