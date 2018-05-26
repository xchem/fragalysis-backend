from graphene_django import DjangoObjectType
import graphene
from graphene_django.rest_framework.mutation import SerializerMutation
from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
)


class Molecule(SerializerMutation):

    class Meta:
        serializer_class = MoleculeSerializer


class Protein(SerializerMutation):

    class Meta:
        serializer_class = ProteinSerializer


class Compound(SerializerMutation):

    class Meta:
        serializer_class = CompoundSerializer


class Target(SerializerMutation):

    class Meta:
        serializer_class = TargetSerializer


class Query(graphene.ObjectType):
    molecules = graphene.List(Molecule)
    proteins = graphene.List(Protein)
    compounds = graphene.List(Compound)
    targets = graphene.List(Target)


schema = graphene.Schema(query=Query)
