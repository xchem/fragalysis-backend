from rest_framework import serializers
from viewer.models import ActivityPoint, Molecule, Project, Protein, Compound, Target, PanddaSite, PanddaEvent, \
    Vector3D,Vector,Interaction,ProteinResidue,TargetResidue, InteractionPoint


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

class Vector3DSerialzier(serializers.ModelSerializer):

    class Meta:
        model = Vector3D
        fields = ('id', 'mol_id', 'vector_id', 'number', 'start_x', 'start_y','start_z', 'end_x', 'end_y', 'end_z')

class VectorSerialzier(serializers.ModelSerializer):

    class Meta:
        model = Vector
        fields = ('id','cmpd_id','smiles','type')

class InteractionSerialzier(serializers.ModelSerializer):

    class Meta:
        model = Interaction
        fields = ('id','interaction_version','interaction_point','interaction_type','distance',
                  'score','prot_smarts','mol_smarts')

class InteractionPointSerializer(serializers.ModelSerializer):

    class Meta:
        model = InteractionPoint
        fields = ('id','prot_res_id', 'mol_id','protein_atom_name','molecule_atom_name',)

class ProteinResidueSerialzier(serializers.ModelSerializer):

    class Meta:
        model = ProteinResidue
        fields = ('id','prot_id','targ_res_id')

class TargetResidueSerialzier(serializers.ModelSerializer):

    class Meta:
        model = TargetResidue
        fields = ('id','target_id','res_name','res_num','chain_id')