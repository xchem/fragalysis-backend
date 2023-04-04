import os
from frag.network.decorate import get_3d_vects_for_mol, get_vect_indices_for_mol
from frag.network.query import get_full_graph
from urllib.parse import urljoin

from api.security import ISpyBSafeQuerySet
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
    ActionType,
    SnapshotActions,
    SessionActions,
    ComputedMolecule,
    ComputedSet,
    NumericalScoreValues,
    ScoreDescription,
    File,
    TextScoreValues,
    TagCategory,
    MoleculeTag,
    SessionProjectTag,
    JobFileTransfer,
    JobRequest
)
from viewer.utils import get_https_host

from scoring.models import MolGroup

from django.contrib.auth.models import User
from django.conf import settings

from rest_framework import serializers

_ISPYB_SAFE_QUERY_SET = ISpyBSafeQuerySet()


class FileSerializer(serializers.ModelSerializer):
    class Meta:
        model = File
        fields = "__all__"

def get_protein_sequences(pdb_block):
    sequence_list = []
    aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
          'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    current_chain = 'A'
    current_sequence = ''
    current_number = 0

    for line in pdb_block.split('\n'):

        if line[0:4] == 'ATOM':
            residue = line[17:20].strip()
            chain = line[21].strip()
            n = int(line[22:26].strip())

            if chain != current_chain:
                chain_dict = {'chain': current_chain, 'sequence': current_sequence}
                sequence_list.append(chain_dict)
                current_sequence = ''
                current_chain = chain
            if not n == current_number:
                if n == current_number + 1:
                    try:
                        seqres = aa[residue]
                    except:
                        seqres = 'X'
                    current_sequence += seqres
                else:
                    if not current_number == 0:
                        gap = n - current_number
                        gap_str = ''
                        for i in range(0, gap):
                            gap_str += 'X'
                        current_sequence += gap_str
            current_number = n

    if not sequence_list:
        chain_dict = {'chain': current_chain, 'sequence': current_sequence}
        sequence_list.append(chain_dict)

    return sequence_list


def protein_sequences(obj):
    """Common enabler code for Target-related serializers
    """
    proteins = obj.protein_set.filter()
    protein_file = None
    for protein in proteins:
        if protein.pdb_info:
            if not os.path.isfile(protein.pdb_info.path):
                continue
            protein_file = protein.pdb_info
            break
    if not protein_file:
        return [{'chain': '', 'sequence': ''}]

    protein_file.open(mode='r')
    pdb_block = protein_file.read()
    protein_file.close()

    sequences = get_protein_sequences(pdb_block)
    return sequences


def template_protein(obj):
    """Common enabler code for Target-related serializers
    """

    proteins = obj.protein_set.filter()
    for protein in proteins:
        if protein.pdb_info:
            return protein.pdb_info.url
    return "NOT AVAILABLE"


class TargetSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()
    zip_archive = serializers.SerializerMethodField()
    metadata = serializers.SerializerMethodField()
    sequences = serializers.SerializerMethodField()

    def get_template_protein(self, obj):
        return template_protein(obj)

    def get_zip_archive(self, obj):
        # The if-check is because the filefield in target has null=True.
        # Note that this link will not work on local
        if hasattr(obj, 'zip_archive') and obj.zip_archive.name:
            return urljoin(get_https_host(self.context["request"]), obj.zip_archive.url)
        return

    def get_metadata(self, obj):
        if hasattr(obj, 'metadata') and obj.metadata.name:
            return urljoin(get_https_host(self.context["request"]), obj.metadata.url)
        return

    def get_sequences(self, obj):
        return protein_sequences(obj)

    class Meta:
        model = Target
        fields = ("id", "title", "project_id", "protein_set", "default_squonk_project",
                  "template_protein", "metadata", "zip_archive", "upload_status", "sequences")


class CompoundSerializer(serializers.ModelSerializer):

    class Meta:
        model = Compound
        fields = (
            "id",
            "inchi",
            "long_inchi",
            "smiles",
            "current_identifier",
            "all_identifiers",
            # "project_id",
            "mol_log_p",
            "mol_wt",
            "tpsa",
            "heavy_atom_count",
            "heavy_atom_mol_wt",
            "nhoh_count",
            "no_count",
            "num_h_acceptors",
            "num_h_donors",
            "num_het_atoms",
            "num_rot_bonds",
            "num_val_electrons",
            "ring_count",
            # "inspirations",
            # "description",
            # "comments",
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
        fields = '__all__'


class ProjectSerializer(serializers.ModelSerializer):

    # Field name translation (prior to refactoring the Model)
    # 'tas' is the new name for 'title'
    target_access_string =  serializers.SerializerMethodField()
    # 'authority' is the (as yet to be implemented) origin of the TAS
    # For now this is fixed at "DIAMOND-ISPYB"
    authority =  serializers.SerializerMethodField()
    # 'can_use_squonk' defines whether a user cna use Squonk for the Projetc
    user_can_use_squonk =  serializers.SerializerMethodField()

    def get_target_access_string(self, instance):
        return instance.title

    def get_user_can_use_squonk(self, instance):
        # User can use Squonk if there is a user object (i.e. they're authenticated)
        # and ISPyB has the user in the Project
        user = self.context['request'].user
        if not user or instance.title not in _ISPYB_SAFE_QUERY_SET.get_proposals_for_user(user):
            return False
        return True

    def get_authority(self, instance):
        # Don't actually need the instance here.
        # We return a hard-coded string.
        del instance
        return "DIAMOND-ISPYB"

    class Meta:
        model = Project
        fields = ("id",
                  "target_access_string",
                  "init_date",
                  "authority",
                  "open_to_public",
                  "user_can_use_squonk")


class MolImageSerializer(serializers.ModelSerializer):

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


class CmpdImageSerializer(serializers.ModelSerializer):

    cmpd_image = serializers.SerializerMethodField()

    def get_cmpd_image(self, obj):
        return draw_mol(obj.smiles, height=125, width=125)

    class Meta:
        model = Compound
        fields = ("id", "cmpd_image")


class ProtMapInfoSerializer(serializers.ModelSerializer):

    map_data = serializers.SerializerMethodField()

    def get_map_data(self, obj):
        if obj.map_info:
            return open(obj.map_info.path, encoding='utf-8').read()
        else:
            return None

    class Meta:
        model = Protein
        fields = ("id", "map_data", "prot_type")


class ProtPDBInfoSerializer(serializers.ModelSerializer):

    pdb_data = serializers.SerializerMethodField()

    def get_pdb_data(self, obj):
        return open(obj.pdb_info.path, encoding='utf-8').read()


    class Meta:
        model = Protein
        fields = ("id", "pdb_data", "prot_type")


class ProtPDBBoundInfoSerializer(serializers.ModelSerializer):

    bound_pdb_data = serializers.SerializerMethodField()

    def get_bound_pdb_data(self, obj):
        if obj.bound_info:
            return open(obj.bound_info.path, encoding='utf-8').read()
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
    graph_auth = os.environ.get("NEO4J_AUTH", "neo4j/neo4j")

    def get_graph(self, obj):
        return get_full_graph(obj.smiles, self.graph_choice, self.graph_auth,
                              isomericSmiles=True)

    class Meta:
        model = Molecule
        fields = ("id", "graph")


### Start of Session Project
class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ('id', 'username', 'email', 'first_name', 'last_name')


# GET
class ActionTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = ActionType
        fields = '__all__'


# GET
class SessionProjectReadSerializer(serializers.ModelSerializer):
    target = TargetSerializer(read_only=True)
    project = ProjectSerializer(read_only=True)
    target = TargetSerializer(read_only=True)
    # This is for the new tags functionality. The old tags field is left here for backwards
    # compatibility
    session_project_tags = serializers.SerializerMethodField()

    def get_session_project_tags(self, obj):
        sp_tags = SessionProjectTag.objects.filter(session_projects=obj.id).values()
        return sp_tags

    class Meta:
        model = SessionProject
        fields = ('id', 'title', 'init_date', 'description',
                  'target', 'project', 'author', 'tags', 'session_project_tags')


# (POST, PUT, PATCH)
class SessionProjectWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = SessionProject
        fields = '__all__'


# (GET, POST, PUT, PATCH)
class SessionActionsSerializer(serializers.ModelSerializer):
    actions = serializers.JSONField()
    class Meta:
        model = SessionActions
        fields = '__all__'


# GET
class SnapshotReadSerializer(serializers.ModelSerializer):
    author  = UserSerializer()
    session_project  = SessionProjectWriteSerializer()
    class Meta:
        model = Snapshot
        fields = ('id', 'type', 'title', 'author', 'description', 'created', 'data',
                  'session_project', 'parent', 'children', 'additional_info')


# (POST, PUT, PATCH)
class SnapshotWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = Snapshot
        fields = ('id', 'type', 'title', 'author', 'description', 'created', 'data',
                  'session_project', 'parent', 'children', 'additional_info')


# (GET, POST, PUT, PATCH)
class SnapshotActionsSerializer(serializers.ModelSerializer):
    actions = serializers.JSONField()
    class Meta:
        model = SnapshotActions
        fields = '__all__'
## End of Session Project


class ComputedSetSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputedSet
        fields = '__all__'


class ComputedMoleculeSerializer(serializers.ModelSerializer):
    # performance issue
    # inspiration_frags = MoleculeSerializer(read_only=True, many=True)
    class Meta:
        model = ComputedMolecule
        fields = '__all__'


class ScoreDescriptionSerializer(serializers.ModelSerializer):
    class Meta:
        model = ScoreDescription
        fields = '__all__'


class NumericalScoreSerializer(serializers.ModelSerializer):
    score = ScoreDescriptionSerializer(read_only=True)
    class Meta:
        model = NumericalScoreValues
        fields = '__all__'


class TextScoreSerializer(serializers.ModelSerializer):
    score = ScoreDescriptionSerializer(read_only=True)
    class Meta:
        model = TextScoreValues
        fields = '__all__'


class ComputedMolAndScoreSerializer(serializers.ModelSerializer):
    numerical_scores = serializers.SerializerMethodField()
    text_scores = serializers.SerializerMethodField()
    pdb_info = serializers.SerializerMethodField()
    # score_descriptions = serializers.SerializerMethodField()

    class Meta:
        model = ComputedMolecule
        fields = (
            "id",
            "sdf_info",
            "name",
            "smiles",
            "pdb_info",
            "compound",
            "computed_set",
            "computed_inspirations",
            "numerical_scores",
            "text_scores",
            # "score_descriptions",
        )

    def get_numerical_scores(self, obj):
        scores = NumericalScoreValues.objects.filter(compound=obj)
        score_dict = {}
        for score in scores:
            score_dict[score.score.name] = score.value
        return score_dict

    def get_text_scores(self, obj):
        scores = TextScoreValues.objects.filter(compound=obj)
        score_dict = {}
        for score in scores:
            score_dict[score.score.name] = score.value
        return score_dict

    def get_pdb_info(self, obj):
        if obj.pdb:
            return obj.pdb.pdb_info.url
        else:
            return None

    # def get_score_descriptions(self, obj):
    #     descriptions = ScoreDescription.objects.filter(computed_set=obj.computed_set)
    #     desc_dict = {}
    #     for desc in descriptions:
    #         desc_dict[desc.name] = desc.description
    #     return desc_dict

# Start of Discourse Serializer

# Class for customer Discourse API
class DiscoursePostWriteSerializer(serializers.Serializer):
    category_name = serializers.CharField(max_length=200)
    parent_category_name = serializers.CharField(max_length=200, initial=settings.DISCOURSE_PARENT_CATEGORY)
    category_colour = serializers.CharField(max_length=10, initial="0088CC")
    category_text_colour = serializers.CharField(max_length=10, initial="FFFFFF")
    post_title = serializers.CharField(max_length=200)
    post_content = serializers.CharField(max_length=2000)
    post_tags = serializers.JSONField()
# End of Discourse Serializer


# Serializer Class for DictToCsv API
class DictToCsvSerializer(serializers.Serializer):
    title = serializers.CharField(max_length=200)
    dict = serializers.DictField()


# Start of Serializers for Tags
class TagCategorySerializer(serializers.ModelSerializer):
    class Meta:
        model = TagCategory
        fields = '__all__'


class MoleculeTagSerializer(serializers.ModelSerializer):
    class Meta:
        model = MoleculeTag
        fields = '__all__'


class SessionProjectTagSerializer(serializers.ModelSerializer):
    class Meta:
        model = SessionProjectTag
        fields = '__all__'


class TargetMoleculesSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()
    zip_archive = serializers.SerializerMethodField()
    metadata = serializers.SerializerMethodField()
    sequences = serializers.SerializerMethodField()
    molecules  = serializers.SerializerMethodField()
    tags_info = serializers.SerializerMethodField()
    tag_categories = serializers.SerializerMethodField()

    def get_template_protein(self, obj):
        return template_protein(obj)

    def get_zip_archive(self, obj):
        # The if-check is because the filefield in target has null=True.
        # Note that this link will not work on local
        if hasattr(obj, 'zip_archive') and obj.zip_archive.name:
            return urljoin(get_https_host(self.context["request"]), obj.zip_archive.url)
        return

    def get_metadata(self, obj):
        if hasattr(obj, 'metadata') and obj.metadata.name:
            return urljoin(get_https_host(self.context["request"]), obj.metadata.url)
        return

    def get_sequences(self, obj):
        return protein_sequences(obj)

    def get_molecules(self, obj):
        mols = Molecule.objects.filter(prot_id__target_id=obj.id)
        molecules = []
        for mol in mols:
            mol_data = {
                'id': mol.id,
                'smiles': mol.smiles,
                'cmpd_id': mol.cmpd_id_id,
                'prot_id': mol.prot_id_id,
                'protein_code': mol.prot_id.code,
                "mol_type": mol.mol_type,
                "molecule_protein": mol.prot_id.pdb_info.url,
                "lig_id": mol.lig_id,
                "chain_id": mol.chain_id,
                'sdf_info': mol.sdf_info,
                'x_com': mol.x_com,
                'y_com': mol.y_com,
                'z_com': mol.z_com,
                'mw': round(mol.cmpd_id.mol_wt, 2),
                'logp': round(mol.cmpd_id.mol_log_p, 2),
                'tpsa': round(mol.cmpd_id.tpsa, 2),
                'ha': mol.cmpd_id.heavy_atom_count,
                'hacc': mol.cmpd_id.num_h_acceptors,
                'hdon': mol.cmpd_id.num_h_donors,
                'rots': mol.cmpd_id.num_rot_bonds,
                'rings': mol.cmpd_id.ring_count,
                'velec': mol.cmpd_id.num_val_electrons
            }
            mol_tags_set = \
                [mt['id'] for mt in MoleculeTag.objects.filter(molecules=mol.id).values()]
            mol_dict = {'data': mol_data, 'tags_set': mol_tags_set}
            molecules.append(mol_dict)

        return molecules

    def get_tags_info(self, obj):
        tags = MoleculeTag.objects.filter(target_id=obj.id)
        tags_info = []
        for tag in tags:
            tag_data = MoleculeTag.objects.filter(id=tag.id).values()
            tag_coords = \
                MolGroup.objects.filter(id=tag.mol_group_id).values('x_com','y_com','z_com' )
            tag_dict = {'data': tag_data, 'coords': tag_coords}
            tags_info.append(tag_dict)

        return tags_info

    def get_tag_categories(self, obj):
        tag_categories = TagCategory.objects.filter(moleculetag__target_id=obj.id).distinct().\
            values()
        return tag_categories

    class Meta:
        model = Target
        fields = ("id", "title", "project_id", "default_squonk_project", "template_protein",
                  "metadata", "zip_archive", "upload_status", "sequences",
                  "molecules", "tags_info", "tag_categories")
# Serializers for Tags - End


# Serializer Class for DownloadStructures API
class DownloadStructuresSerializer(serializers.Serializer):
    target_name = serializers.CharField(max_length=200)
    proteins = serializers.CharField(max_length=5000)
    pdb_info = serializers.BooleanField(default=False)
    bound_info = serializers.BooleanField(default=False)
    cif_info = serializers.BooleanField(default=False)
    mtz_info = serializers.BooleanField(default=False)
    diff_info = serializers.BooleanField(default=False)
    event_info = serializers.BooleanField(default=False)
    sigmaa_info = serializers.BooleanField(default=False)
    sdf_info = serializers.BooleanField(default=False)
    single_sdf_file = serializers.BooleanField(default=False)
    trans_matrix_info = serializers.BooleanField(default=False)
    metadata_info = serializers.BooleanField(default=False)
    smiles_info = serializers.BooleanField(default=False)
    static_link = serializers.BooleanField(default=False)
    file_url = serializers.CharField(max_length=200)


# Start of Serializers for Squonk Jobs
# (GET)
class JobFileTransferReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobFileTransfer
        fields = '__all__'


# (POST, PUT, PATCH)
class JobFileTransferWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobFileTransfer
        fields = ("snapshot", "target", "squonk_project",
                  "proteins", "compounds")


# (GET)
class JobRequestReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobRequest
        fields = '__all__'


# (POST, PUT, PATCH)
class JobRequestWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobRequest
        fields = ("squonk_job_name", "snapshot", "target", "squonk_project",
                  "squonk_job_spec")


# Contains output fields from Fragalysis
class JobCallBackReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = JobRequest
        fields = ("job_status", "job_status_datetime",
                  "squonk_job_name", "squonk_job_spec", "upload_status")


# Contains the input fields from Squonk
class JobCallBackWriteSerializer(serializers.ModelSerializer):
    SQUONK_STATUS = ['PENDING', 'STARTED', 'SUCCESS', 'FAILURE', 'RETRY', 'REVOKED']
    job_status = serializers.ChoiceField(choices=SQUONK_STATUS, default="PENDING")
    state_transition_time = serializers.DateTimeField(source='job_status_datetime')

    class Meta:
        model = JobRequest
        fields = ("job_status", "state_transition_time")

# End of Serializers for Squonk Jobs
