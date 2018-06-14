from rest_framework import serializers
from viewer.models import ActivityPoint, Molecule, Project, Protein, Compound, Target


class TargetSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()
    def get_template_protein(self, obj):
        return obj.protein_set.filter()[0].pdb_info.url


    class Meta:
        model = Target
        fields = ("id", "title", "project_id", "protein_set","template_protein")


class CompoundSerializer(serializers.ModelSerializer):

    class Meta:
        model = Compound
        fields = (
            "id",
            "inchi",
            "smiles",
            "mol_log_p",
            "mol_wt",
            "num_h_acceptors",
            "num_h_donors",
        )


class MoleculeSerializer(serializers.ModelSerializer):

    molecule_protein = serializers.SerializerMethodField()
    def get_molecule_protein(self, obj):
        return obj.prot_id.pdb_info.url

    class Meta:
        model = Molecule
        fields = (
            "id",
            "smiles",
            "cmpd_id",
            "prot_id",
            "molecule_protein"
            "lig_id",
            "chain_id",
            "sdf_info",
            "x_com",
            "y_com",
            "z_com",
        )


class ActivityPointSerializer(serializers.ModelSerializer):

    class Meta:
        model = ActivityPoint
        fields = (
            "id",
            "source",
            "target_id",
            "cmpd_id",
            "activity",
            "units",
            "confidence",
            "operator",
            "internal_id",
        )


class ProteinSerializer(serializers.ModelSerializer):

    class Meta:
        model = Protein
        fields = (
            "id",
            "code",
            "target_id",
            "pdb_info",
            "mtz_info",
            "map_info",
            "cif_info",
        )


class ProjectSerializer(serializers.ModelSerializer):

    class Meta:
        model = Project
        fields = ("id", "title")
