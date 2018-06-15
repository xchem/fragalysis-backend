from rest_framework import serializers
from viewer.models import ActivityPoint, Molecule, Project, Protein, Compound, Target
from api.utils import draw_mol
from frag.network.decorate import get_3d_vects_for_mol
from frag.network.query import get_full_graph


class TargetSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()

    def get_template_protein(self, obj):
        if len(obj.protein_set.filter()) > 0:
            return obj.protein_set.filter()[0].pdb_info.url
        else:
            return "NOT AVAILABLE"

    class Meta:
        model = Target
        fields = ("id", "title", "project_id", "protein_set", "template_protein")


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
            "molecule_protein",
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


class MolImageSerialzier(serializers.ModelSerializer):

    mol_image = serializers.SerializerMethodField()

    def get_mol_image(self, obj):
        return draw_mol(obj.smiles, height=125, width=125)

    class Meta:
        model = Molecule
        fields = ("id", "mol_image")


class CmpdImageSerialzier(serializers.ModelSerializer):

    cmpd_image = serializers.SerializerMethodField()

    def get_cmpd_image(self, obj):
        return draw_mol(obj.smiles, height=125, width=125)

    class Meta:
        model = Compound
        fields = ("id", "cmpd_image")


class ProtMapInfoSerialzer(serializers.ModelSerializer):

    map_data = serializers.SerializerMethodField()

    def get_map_data(self, obj):
        return obj.map_info

    class Meta:
        model = Protein
        fields = ("id", "map_data")


class ProtPDBInfoSerialzer(serializers.ModelSerializer):

    pdb_data = serializers.SerializerMethodField()

    def get_pdb_data(self, obj):
        return open(obj.pdb_info.path).read()

    class Meta:
        model = Protein
        fields = ("id", "pdb_data")


class VectorsSerializer(serializers.ModelSerializer):

    vectors = serializers.SerializerMethodField()

    def get_vectors(self, obj):
        return get_3d_vects_for_mol(obj.sdf_info)

    class Meta:
        model = Protein
        fields = ("id", "vectors")


class GraphSerializer(serializers.ModelSerializer):

    graph = serializers.SerializerMethodField()

    def get_graph(self, obj):
        return get_full_graph(obj.smiles)

    class Meta:
        model = Protein
        fields = ("id", "graph")
