from viewer.models import ViewScene, ActivityPoint, Molecule, Project, Protein, Compound,Target
from viewer.serializers import ViewSceneSerializer,ActivityPointSerializer,MoleculeSerializer,ProjectSerializer,\
    ProteinSerializer, CompoundSerializer, TargetSerializer
from rest_framework import generics
from rest_framework import permissions
from rest_framework import viewsets

class TargetView(generics.ListAPIView,
                 viewsets.ModelViewSet):
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    # permission_classes =  [permissions.DjangoObjectPermissions(),]
    filter_fields = ('title',)

class MoleculeView(generics.ListAPIView):
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    filter_fields = ('prot_id', 'cmpd_id','smiles')

class CompoundView(generics.ListAPIView,
                   viewsets.ModelViewSet):
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_fields = ('smiles',)

class ProteinView(generics.ListAPIView,
                  viewsets.ModelViewSet):
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    filter_fields = ('code','target_id',)
