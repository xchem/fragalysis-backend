from graphene_django import DjangoObjectType
import graphene
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


class Vector(SerializerMutation):

    class Meta:
        serializer_class = VectorSerializer


class InteractionPoint(SerializerMutation):

    class Meta:
        serializer_class = InteractionPointSerializer


class Interaction(SerializerMutation):

    class Meta:
        serializer_class = InteractionSerializer


class ProteinResidue(SerializerMutation):

    class Meta:
        serializer_class = ProteinResidueSerialzier


class TargetResidue(SerializerMutation):

    class Meta:
        serializer_class = TargetResidueSerialzier


class Query(graphene.ObjectType):
    vector3d = graphene.List(Vector3D)
    vector = graphene.List(Vector)
    interaction_point = graphene.List(InteractionPoint)
    interaction = graphene.List(Interaction)
    protein_residue = graphene.List(ProteinResidue)
    target_residue = graphene.List(TargetResidue)


schema = graphene.Schema(query=Query)
