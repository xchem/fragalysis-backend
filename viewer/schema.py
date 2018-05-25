from graphene_django import DjangoObjectType
import graphene
from graphene_django.rest_framework.mutation import SerializerMutation
from viewer.models import Molecule, Project, Protein, Compound, Target
from viewer.serializers import (
    MoleculeSerializer,
    ProteinSerializer,
    CompoundSerializer,
    TargetSerializer,
)


class Molecule(SerializerMutation):

    class Meta:
        serializer_class = MoleculeSerializer


class Query(graphene.ObjectType):
    molecules = graphene.List(Molecule)

    def resolve_users(self, info):
        return Molecule.objects.all()


schema = graphene.Schema(query=Query)
