import graphene
from graphene import relay
from graphene_django.rest_framework.mutation import SerializerMutation

from hypothesis.serializers import (
    Vector3DSerializer,
    VectorSerializer,
    InteractionPointSerializer,
    InteractionSerializer,
    ProteinResidueSerialzier,
    TargetResidueSerialzier,
)


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


class ProteinResidue(SerializerMutation):

    class Meta:
        serializer_class = ProteinResidueSerialzier
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
    protein_residue = graphene.List(ProteinResidue)
    target_residue = graphene.List(TargetResidue)


schema = graphene.Schema(query=Query)
