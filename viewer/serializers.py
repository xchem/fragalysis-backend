from rest_framework import serializers
from viewer.models import ViewScene, ActivityPoint, Molecule, Project, Protein, Compound, Target


class TargetSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Target
        fields = ('title', 'project_id')

class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = Compound
        fields = ('inchi', 'smiles', 'mol_log_p', 'mol_wt', 'num_h_acceptors', 'num_h_donors',)

class MoleculeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Molecule
        fields = ('smiles', 'cmpd_id', 'prot_id', 'lig_id', 'chain_id', 'sdf_info', 'x_com', 'y_com', 'z_com',)

class ActivityPointSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = ActivityPoint
        fields = ('source','target_id', 'cmpd_id', 'activity', 'units', 'confidence', 'operator', 'internal_id')

class ProteinSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Protein
        fields = ('code','target_id','pdb_info','mtz_info','map_info','cif_info')

class ProjectSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Project
        fields = ('title',)

class ViewSceneSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model =  ViewScene
        fields = ('uuid','title','scene',)