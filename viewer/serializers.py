from rest_framework import serializers
from .models import ViewScene, ActivityPoint, Molecule, Project, Protein, Compound

class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = Compound
        fields = ('inchi', 'smiles', 'mol_log_p', 'mol_wt', 'num_h_acceptors', 'num_h_donors',)

class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Molecule
        fields = ('smiles', 'cmpd_id', 'prot_id', 'lig_id', 'chain_id', 'sdf_info', 'x_com', 'y_com', 'z_com',)

class ActivityPointSerializer(serializers.HyperlinkedIdentityField):
    class Meta:
        model = ActivityPoint
        fields = ('source','target_id', 'cmpd_id', 'activity', 'units', 'confidence', 'operator', 'internal_id')

class ProteinSerializer(serializers.HyperlinkedIdentityField):
    class Meta:
        model = Protein
        fields = ('code','target_id','pdb_info','cif_info','mtz_info',)

class ProjectSerializer(serializers.HyperlinkedIdentityField):
    class Meta:
        model = Project
        fields = ('title',)

class ViewSceneSerializer(serializers.HyperlinkedIdentityField):
    class Meta:
        model =  ViewScene
        fields = ('uuid','title','scene',)