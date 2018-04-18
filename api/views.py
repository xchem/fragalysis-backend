from viewer.models import ActivityPoint, Molecule, Project, Protein, Compound,Target, PanddaEvent, PanddaSite
from viewer.serializers import MoleculeSerializer, ProteinSerializer, CompoundSerializer, TargetSerializer, \
    PanddaEventSerializer, PanddaSiteSerializer
from rest_framework import permissions
from rest_framework import viewsets


class TargetView(viewsets.ReadOnlyModelViewSet):
    queryset = Target.objects.filter()
    serializer_class = TargetSerializer
    # permission_classes =  [permissions.DjangoObjectPermissions,]
    filter_fields = ('title',)

class MoleculeView(viewsets.ReadOnlyModelViewSet):
    queryset = Molecule.objects.filter()
    serializer_class = MoleculeSerializer
    filter_fields = ('prot_id', 'cmpd_id','smiles','prot_id__target_id', 'mol_groups')

class CompoundView(viewsets.ReadOnlyModelViewSet):
    queryset = Compound.objects.filter()
    serializer_class = CompoundSerializer
    filter_fields = ('smiles',)

class ProteinView(viewsets.ReadOnlyModelViewSet):
    queryset = Protein.objects.filter()
    serializer_class = ProteinSerializer
    filter_fields = ('code','target_id',)

class PanddaEventView(viewsets.ReadOnlyModelViewSet):
    queryset = PanddaEvent.objects.filter()
    serializer_class = PanddaEventSerializer
    filter_fields = ('xtal','event','pandda_site','target_id',)

class PanddaSiteView(viewsets.ReadOnlyModelViewSet):
    queryset = PanddaSite.objects.filter()
    serializer_class = PanddaSiteSerializer
    filter_fields = ('pandda_run','site_id','pandda_version','target_id',)

# TODO VIEW to get Images, PDBs, Vectors and Graphs

"""
    # New REST functions
    url(r'^img_from_mol_pk/(?P<pk>[0-9]+)/$', views.img_from_mol_pk, name='img_from_mol_pk'),
    url(r'^img_from_cmpd_pk/(?P<pk>[0-9]+)/$', views.img_from_cmpd_pk, name='img_from_cmpd_pk'),

    # Move all of these to API
    url(r'^get_vects_from_pk/(?P<pk>[0-9]+)/', views.get_vects_from_pk, name='get_vects_from_pk'),
    url(r'^get_graph_from_pk/(?P<pk>[0-9]+)/', views.get_graph_from_pk, name='get_graph_from_pk'),
"""
