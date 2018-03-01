from viewer.models import ViewScene, ActivityPoint, Molecule, Project, Protein, Compound,Target
from viewer.serializers import ViewSceneSerializer,ActivityPointSerializer,MoleculeSerializer,ProjectSerializer,\
    ProteinSerializer, CompoundSerializer, TargetSerializer, MDLSerializer
from rest_framework import permissions
from rest_framework import viewsets


class TargetView(viewsets.ModelViewSet):
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    # permission_classes =  [permissions.DjangoObjectPermissions,]
    filter_fields = ('title',)

class MoleculeView(viewsets.ModelViewSet):
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    filter_fields = ('prot_id', 'cmpd_id','smiles','prot_id__target_id')

class MDLView(viewsets.ModelViewSet):
    queryset = Molecule.objects.filter()
    serializer_class = MDLSerializer

class CompoundView(viewsets.ModelViewSet):
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_fields = ('smiles',)

class ProteinView(viewsets.ModelViewSet):
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    filter_fields = ('code','target_id',)
