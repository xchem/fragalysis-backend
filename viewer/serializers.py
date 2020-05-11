import os
from frag.network.decorate import get_3d_vects_for_mol, get_vect_indices_for_mol
from frag.network.query import get_full_graph
from rest_framework import serializers

from api.utils import draw_mol
from viewer.models import (
    ActivityPoint,
    Molecule,
    Project,
    Protein,
    Compound,
    Target,
    Snapshot,
    SessionProject,
    ComputedCompound,
    CompoundSet,
    NumericalScoreValues,
    ScoreDescription
)
from django.contrib.auth.models import User

class TargetSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()

    def get_template_protein(self, obj):
        proteins = obj.protein_set.filter()
        for protein in proteins:
            if protein.pdb_info:
                return protein.pdb_info.url
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
    protein_code = serializers.SerializerMethodField()
    mw = serializers.SerializerMethodField()
    logp = serializers.SerializerMethodField()
    tpsa = serializers.SerializerMethodField()
    ha = serializers.SerializerMethodField()
    hacc = serializers.SerializerMethodField()
    hdon = serializers.SerializerMethodField()
    rots = serializers.SerializerMethodField()
    rings = serializers.SerializerMethodField()
    velec = serializers.SerializerMethodField()

    def get_molecule_protein(self, obj):
        return obj.prot_id.pdb_info.url

    def get_protein_code(self, obj):
        return obj.prot_id.code

    def get_mw(self, obj):
        return round(obj.cmpd_id.mol_wt, 2)

    def get_logp(self, obj):
        return round(obj.cmpd_id.mol_log_p, 2)

    def get_tpsa(self, obj):
        return round(obj.cmpd_id.tpsa, 2)

    def get_ha(self, obj):
        return obj.cmpd_id.heavy_atom_count

    def get_hacc(self, obj):
        return obj.cmpd_id.num_h_acceptors

    def get_hdon(self, obj):
        return obj.cmpd_id.num_h_donors

    def get_rots(self, obj):
        return obj.cmpd_id.num_rot_bonds

    def get_rings(self, obj):
        return obj.cmpd_id.ring_count

    def get_velec(self, obj):
        return obj.cmpd_id.num_val_electrons

    class Meta:
        model = Molecule
        fields = (
            "id",
            "smiles",
            "cmpd_id",
            "prot_id",
            "protein_code",
            "mol_type",
            "molecule_protein",
            "lig_id",
            "chain_id",
            "sdf_info",
            "x_com",
            "y_com",
            "z_com",
            "mw",
            "logp",
            "tpsa",
            "ha",
            "hacc",
            "hdon",
            "rots",
            "rings",
            "velec"
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
            "prot_type",
            "pdb_info",
            "bound_info",
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
        request = self.context["request"]
        params = request.query_params
        if params:
            return draw_mol(
                obj.smiles,
                height=int(float(params["height"])),
                width=int(float(params["width"])),
            )
        else:
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
        if obj.map_info:
            return open(obj.map_info.path).read()
        else:
            return None

    class Meta:
        model = Protein
        fields = ("id", "map_data", "prot_type")


class ProtPDBInfoSerialzer(serializers.ModelSerializer):

    pdb_data = serializers.SerializerMethodField()

    def get_pdb_data(self, obj):
        return open(obj.pdb_info.path).read()


    class Meta:
        model = Protein
        fields = ("id", "pdb_data", "prot_type")


class ProtPDBBoundInfoSerialzer(serializers.ModelSerializer):

    bound_pdb_data = serializers.SerializerMethodField()

    def get_bound_pdb_data(self, obj):
        if obj.bound_info:
            return open(obj.bound_info.path).read()
        else:
            return None

    class Meta:
        model = Protein
        fields = ("id", "bound_pdb_data", "target_id")


class VectorsSerializer(serializers.ModelSerializer):

    vectors = serializers.SerializerMethodField()

    def get_vectors(self, obj):
        out_data = {}
        try:
            out_data["3d"] = get_3d_vects_for_mol(obj.sdf_info, iso_labels=False)
        # temporary patch
        except IndexError:
            out_data["3d"] = get_3d_vects_for_mol(obj.sdf_info, iso_labels=True)
        out_data["indices"] = get_vect_indices_for_mol(obj.sdf_info)
        return out_data

    class Meta:
        model = Molecule
        fields = ("id", "vectors")


class GraphSerializer(serializers.ModelSerializer):

    graph = serializers.SerializerMethodField()
    graph_choice = os.environ.get("NEO4J_QUERY", "neo4j")

    def get_graph(self, obj):
        return get_full_graph(obj.smiles, self.graph_choice, isomericSmiles=True)

    class Meta:
        model = Molecule
        fields = ("id", "graph")


### Start of Session Project
class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ('id', 'username', 'email', 'first_name', 'last_name')


# GET
class SessionProjectReadSerializer(serializers.ModelSerializer):
    target = TargetSerializer(read_only=True)
    author = UserSerializer(read_only=True)
    class Meta:
        model = SessionProject
        fields = '__all__'

# (POST, PUT, PATCH)
class SessionProjectWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = SessionProject
        fields = '__all__'


# GET
class SnapshotReadSerializer(serializers.ModelSerializer):
    author  = UserSerializer()
    session_project  = SessionProjectWriteSerializer()
    class Meta:
        model = Snapshot
        fields = ('id', 'type', 'title', 'author', 'description', 'created', 'data', 'session_project', 'parent', 'children')


# (POST, PUT, PATCH)
class SnapshotWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = Snapshot
        fields = ('id', 'type', 'title', 'author', 'description', 'created', 'data', 'session_project', 'parent', 'children')
## End of Session Project

# class FileSerializer(serializers.ModelSerializer):
#     class Meta:
#         model = ComputedCompound
#         fields = "__all__"

class ScoreDescriptionSerializer(serializers.Serializer):
    name = serializers.CharField()
    description = serializers.CharField()

class TargetCompoundSetsSerializer(serializers.Serializer):
    score = ScoreDescriptionSerializer()
    value = serializers.FloatField()


class CompoundSetSerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundSet
        fields = '__all__'


class CompoundMoleculeSerializer(serializers.ModelSerializer):
    # performance issue
    # inspiration_frags = MoleculeSerializer(read_only=True, many=True)
    class Meta:
        model = ComputedCompound
        fields = '__all__'



class ScoreDescriptionAllSerializer(serializers.ModelSerializer):
    class Meta:
        model = ScoreDescription
        fields = '__all__'

class NumericalScoreSerializer(serializers.ModelSerializer):
    score = ScoreDescriptionAllSerializer(read_only=True)
    class Meta:
        model = NumericalScoreValues
        fields = '__all__'