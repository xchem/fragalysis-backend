import graphene
from graphene_django.rest_framework.mutation import SerializerMutation
from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
)

relay = graphene.relay


class Molecule(SerializerMutation):

    class Meta:
        serializer_class = MoleculeSerializer
        interfaces = (relay.Node,)


class Protein(SerializerMutation):

    class Meta:
        serializer_class = ProteinSerializer
        interfaces = (relay.Node,)


class Compound(SerializerMutation):

    class Meta:
        serializer_class = CompoundSerializer
        interfaces = (relay.Node,)


class Target(SerializerMutation):

    class Meta:
        serializer_class = TargetSerializer
        interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    molecules = graphene.List(Molecule)
    proteins = graphene.List(Protein)
    compounds = graphene.List(Compound)
    targets = graphene.List(Target)


schema = graphene.Schema(query=Query)
