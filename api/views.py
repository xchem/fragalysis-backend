from viewer.models import ViewScene, ActivityPoint, Molecule, Project, Protein, Compound,Target
from viewer.serializers import ViewSceneSerializer,ActivityPointSerializer,MoleculeSerializer,ProjectSerializer,\
    ProteinSerializer, CompoundSerializer
from rest_framework import generics
from rest_framework import permissions

class TargetList(generics.ListAPIView):
    queryset = Target.objects.filter()
    serializer_class = CompoundSerializer
    permission_classes =  [permissions.DjangoObjectPermissions(),]
    filter_fields = ('title',)

class MoleculeList(generics.ListAPIView):
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    permission_classes =  [permissions.DjangoObjectPermissions(),]
    filter_fields = ('prot_id', 'cmpd_id','smiles')

class CompoundList(generics.ListAPIView):
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    permission_classes =  [permissions.DjangoObjectPermissions(),]
    filter_fields = ('smiles',)

class ProteinList(generics.ListAPIView):
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    permission_classes =  [permissions.DjangoObjectPermissions(),]
    filter_fields = ('code','target_id',)
