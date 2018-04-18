from rest_framework import serializers
from viewer.models import ActivityPoint, Molecule, Project, Protein, Compound, Target, PanddaSite, PanddaEvent


class TargetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Target
        fields = ('id', 'title', 'project_id','protein_set')

class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = Compound
        fields = ('id', 'inchi', 'smiles', 'mol_log_p', 'mol_wt', 'num_h_acceptors', 'num_h_donors',)

class MoleculeSerializer(serializers.ModelSerializer):
    class Meta:
        model = Molecule
        fields = ('id','smiles', 'cmpd_id', 'prot_id', 'lig_id', 'chain_id', 'sdf_info', 'x_com', 'y_com', 'z_com',)

class ActivityPointSerializer(serializers.ModelSerializer):
    class Meta:
        model = ActivityPoint
        fields = ('id', 'source','target_id', 'cmpd_id', 'activity', 'units', 'confidence', 'operator', 'internal_id')

class ProteinSerializer(serializers.ModelSerializer):

    class Meta:
        model = Protein
        fields = ('id', 'code','target_id','pdb_info','mtz_info','map_info','cif_info')

class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = ('id', 'title',)

class PanddaEventSerializer(serializers.ModelSerializer):

    class Meta:
        model = PanddaEvent
        fields = ('id','xtal','event','pandda_site','target_id','pdb_info','mtz_info','map_info','lig_id',
                  'event_com_x','event_com_y','event_com_z','lig_com_x','lig_com_y','lig_com_z',
                  'event_dist_from_site_centroid','lig_dist_from_site_centroid','small_map_info')

class PanddaSiteSerializer(serializers.ModelSerializer):

    class Meta:
        model = PanddaSite
        fields = ('id','site_id','pandda_run','pandda_version','target_id','site_align_com_x','site_align_com_y','site_align_com_z',
                  'site_native_com_x','site_native_com_y','site_native_com_z',)