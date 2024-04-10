import logging

# import os
from pathlib import Path
from urllib.parse import urljoin

import yaml
from django.conf import settings
from django.contrib.auth.models import User

# from django.contrib.postgres.aggregates import ArrayAgg
# from django.db.models import F
from frag.network.decorate import get_3d_vects_for_mol, get_vect_indices_for_mol
from frag.network.query import get_full_graph
from rdkit import Chem
from rdkit.Chem import Descriptors
from rest_framework import serializers

from api.security import ISpyBSafeQuerySet
from api.utils import draw_mol, validate_tas

# from scoring.models import SiteObservationGroup
from viewer import models
from viewer.target_loader import XTALFORMS_FILE
from viewer.target_set_upload import sanitize_mol
from viewer.utils import get_https_host

logger = logging.getLogger(__name__)

_ISPYB_SAFE_QUERY_SET = ISpyBSafeQuerySet()


class FileSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.File
        fields = "__all__"


# def get_protein_sequences(pdb_block):
#     sequence_list = []
#     aa = {
#         'CYS': 'C',
#         'ASP': 'D',
#         'SER': 'S',
#         'GLN': 'Q',
#         'LYS': 'K',
#         'ILE': 'I',
#         'PRO': 'P',
#         'THR': 'T',
#         'PHE': 'F',
#         'ASN': 'N',
#         'GLY': 'G',
#         'HIS': 'H',
#         'LEU': 'L',
#         'ARG': 'R',
#         'TRP': 'W',
#         'ALA': 'A',
#         'VAL': 'V',
#         'GLU': 'E',
#         'TYR': 'Y',
#         'MET': 'M',
#     }

#     current_chain = 'A'
#     current_sequence = ''
#     current_number = 0

#     for line in pdb_block.split('\n'):
#         if line[0:4] == 'ATOM':
#             residue = line[17:20].strip()
#             chain = line[21].strip()
#             n = int(line[22:26].strip())

#             if chain != current_chain:
#                 chain_dict = {'chain': current_chain, 'sequence': current_sequence}
#                 sequence_list.append(chain_dict)
#                 current_sequence = ''
#                 current_chain = chain
#             if not n == current_number:
#                 if n == current_number + 1:
#                     try:
#                         seqres = aa[residue]
#                     except:
#                         seqres = 'X'
#                     current_sequence += seqres
#                 else:
#                     if not current_number == 0:
#                         gap = n - current_number
#                         gap_str = ''
#                         for _ in range(gap):
#                             gap_str += 'X'
#                         current_sequence += gap_str
#             current_number = n

#     if not sequence_list:
#         chain_dict = {'chain': current_chain, 'sequence': current_sequence}
#         sequence_list.append(chain_dict)

#     return sequence_list


# def protein_sequences(obj):
#     """Common enabler code for Target-related serializers"""
#     proteins = models.SiteObservation.filter_manager.by_target(target=obj)
#     protein_file = None
#     for protein in proteins:
#         if protein.apo_file:
#             if not os.path.isfile(protein.apo_file.path):
#                 continue
#             protein_file = protein.apo_file
#             break
#     if not protein_file:
#         return [{'chain': '', 'sequence': ''}]

#     protein_file.open(mode='r')
#     pdb_block = protein_file.read()
#     protein_file.close()

#     sequences = get_protein_sequences(pdb_block)
#     return sequences


# def template_protein(obj):
#     """Common enabler code for Target-related serializers"""

#     proteins = models.SiteObservation.filter_manager.by_target(target=obj)
#     for protein in proteins:
#         if protein.apo_file:
#             return protein.apo_file.url
#     return "NOT AVAILABLE"


class CompoundIdentifierTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CompoundIdentifierType
        fields = '__all__'


class CompoundIdentifierSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CompoundIdentifier
        fields = '__all__'


class TargetSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()
    zip_archive = serializers.SerializerMethodField()
    metadata = serializers.SerializerMethodField()
    # sequences = serializers.SerializerMethodField()

    def get_template_protein(self, obj):
        exp_upload = (
            models.ExperimentUpload.objects.filter(
                target=obj,
            )
            .order_by('-commit_datetime')
            .first()
        )

        yaml_path = (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(exp_upload.task_id)
        )

        # add unpacked zip directory
        yaml_path = [d for d in list(yaml_path.glob("*")) if d.is_dir()][0]

        # add upload_[d] dir
        yaml_path = next(yaml_path.glob("upload_*"))

        # last components of path, need for reconstruction later
        comps = yaml_path.parts[-2:]

        # and the file itself
        yaml_path = yaml_path.joinpath(XTALFORMS_FILE)
        logger.debug("assemblies path: %s", yaml_path)
        if yaml_path.is_file():
            with open(yaml_path, "r", encoding="utf-8") as file:
                contents = yaml.safe_load(file)
                try:
                    assemblies = contents["assemblies"]
                except KeyError:
                    logger.error("No 'assemblies' section in '%s'", XTALFORMS_FILE)
                    return ''

                try:
                    first = list(assemblies.values())[0]
                except IndexError:
                    logger.error("No assemblies in 'assemblies' section")
                    return ''

                try:
                    reference = first["reference"]
                except KeyError:
                    logger.error("No assemblies in 'assemblies' section")
                    return ''

                ref_path = (
                    Path(settings.TARGET_LOADER_MEDIA_DIRECTORY)
                    .joinpath(exp_upload.task_id)
                    .joinpath(comps[0])
                    .joinpath(comps[1])
                    .joinpath("crystallographic_files")
                    .joinpath(reference)
                    .joinpath(f"{reference}.pdb")
                )
                logger.debug('ref_path: %s', ref_path)
                if Path(settings.MEDIA_ROOT).joinpath(ref_path).is_file():
                    request = self.context.get('request', None)
                    if request is not None:
                        return request.build_absolute_uri(
                            Path(settings.MEDIA_URL).joinpath(ref_path)
                        )
                    else:
                        return ''
                else:
                    logger.error("Reference pdb file doesn't exist")
                    return ''
        else:
            logger.error("'%s' missing", XTALFORMS_FILE)
            return ''

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

    # def get_sequences(self, obj):
    #     return protein_sequences(obj)

    class Meta:
        model = models.Target
        # TODO: it's missing protein_set. is it necessary anymore?
        fields = (
            "id",
            "title",
            "display_name",
            "project_id",
            "default_squonk_project",
            "template_protein",
            "metadata",
            "zip_archive",
            "upload_status",
            # "sequences",
        )
        extra_kwargs = {
            "id": {"read_only": True},
            "title": {"read_only": True},
            "project_id": {"read_only": True},
            "default_squonk_project": {"read_only": True},
            "template_protein": {"read_only": True},
            "metadata": {"read_only": True},
            "zip_archive": {"read_only": True},
            "upload_status": {"read_only": True},
            # "sequences": {"read_only": True},
        }


class CompoundSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Compound
        fields = (
            "id",
            "inchi",
            "smiles",
            "current_identifier",
            "all_identifiers",
            "project_id",
            "compound_code",
            "inspirations",
            "description",
            "comments",
        )


class SiteObservationSerializer(serializers.ModelSerializer):
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
        return obj.experiment.pdb_info.url

    def get_mw(self, obj):
        return round(obj.cmpd.mol_wt, 2)

    def get_logp(self, obj):
        return round(obj.cmpd.mol_log_p, 2)

    def get_tpsa(self, obj):
        return round(obj.cmpd.tpsa, 2)

    def get_ha(self, obj):
        return obj.cmpd.heavy_atom_count

    def get_hacc(self, obj):
        return obj.cmpd.num_h_acceptors

    def get_hdon(self, obj):
        return obj.cmpd.num_h_donors

    def get_rots(self, obj):
        return obj.cmpd.num_rot_bonds

    def get_rings(self, obj):
        return obj.cmpd.ring_count

    def get_velec(self, obj):
        return obj.cmpd.num_val_electrons

    class Meta:
        model = models.SiteObservation
        fields = (
            "id",
            "smiles",
            "cmpd",
            "code",
            # seems that this field is removed
            # "mol_type",
            "molecule_protein",
            "seq_id",
            "chain_id",
            "ligand_mol_file",
            # these seem to be removed as well
            # "x_com",
            # "y_com",
            # "z_com",
            "mw",
            "logp",
            "tpsa",
            "ha",
            "hacc",
            "hdon",
            "rots",
            "rings",
            "velec",
        )


class ActivityPointSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ActivityPoint
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


class ProjectSerializer(serializers.ModelSerializer):
    # Field name translation (prior to refactoring the Model)
    # 'tas' is the new name for 'title'
    target_access_string = serializers.SerializerMethodField()
    # 'authority' is the (as yet to be implemented) origin of the TAS
    # For now this is fixed at "DIAMOND-ISPYB"
    authority = serializers.SerializerMethodField()
    # 'can_use_squonk' defines whether a user cna use Squonk for the Projetc
    user_can_use_squonk = serializers.SerializerMethodField()

    def get_target_access_string(self, instance):
        return instance.title

    def get_user_can_use_squonk(self, instance):
        # User can use Squonk if there is a user object (i.e. they're authenticated)
        # and ISPyB has the user in the Project
        user = self.context['request'].user
        if (
            not user
            or instance.title not in _ISPYB_SAFE_QUERY_SET.get_proposals_for_user(user)
        ):
            return False
        return True

    def get_authority(self, instance):
        # Don't actually need the instance here.
        # We return a hard-coded string.
        del instance
        return "DIAMOND-ISPYB"

    class Meta:
        model = models.Project
        fields = (
            "id",
            "target_access_string",
            "init_date",
            "authority",
            "open_to_public",
            "user_can_use_squonk",
        )


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
        model = models.SiteObservation
        fields = ("id", "mol_image")


class CmpdImageSerializer(serializers.ModelSerializer):
    cmpd_image = serializers.SerializerMethodField()

    def get_cmpd_image(self, obj):
        return draw_mol(obj.smiles, height=125, width=125)

    class Meta:
        model = models.Compound
        fields = ("id", "cmpd_image")


# TODO: this is a tricky one. there's event_map_info in Experiment,
# but this is array_field of files. there's no way to link back to a
# single protein, I can only do target now
class ProtMapInfoSerializer(serializers.ModelSerializer):
    map_data = serializers.SerializerMethodField()

    def get_map_data(self, obj):
        # TODO: this was formerly protein.map_info. based on
        # description in the spec, it seems this is equivalent to
        # event_file in site_observation, but it's binary and not read
        if obj.event_file:
            return open(obj.event_file.path, encoding='utf-8').read()
        else:
            return None

    class Meta:
        model = models.SiteObservation
        fields = ("id", "map_data")


class ProtPDBInfoSerializer(serializers.ModelSerializer):
    pdb_data = serializers.SerializerMethodField()

    def get_pdb_data(self, obj):
        return open(obj.experiment.pdb_info.path, encoding='utf-8').read()

    class Meta:
        model = models.SiteObservation
        fields = ("id", "pdb_data")


class ProtPDBBoundInfoSerializer(serializers.ModelSerializer):
    bound_pdb_data = serializers.SerializerMethodField()
    target = serializers.IntegerField()

    def get_bound_pdb_data(self, obj):
        if obj.bound_file:
            return open(obj.bound_file.path, encoding='utf-8').read()
        else:
            return None

    class Meta:
        model = models.SiteObservation
        fields = ("id", "bound_pdb_data", "target")


class VectorsSerializer(serializers.ModelSerializer):
    vectors = serializers.SerializerMethodField()

    def get_vectors(self, obj):
        out_data = {}
        if obj.ligand_mol_file:
            try:
                out_data["3d"] = get_3d_vects_for_mol(
                    obj.ligand_mol_file, iso_labels=False
                )
            # temporary patch
            except IndexError:
                out_data["3d"] = get_3d_vects_for_mol(
                    obj.ligand_mol_file, iso_labels=True
                )
            out_data["indices"] = get_vect_indices_for_mol(obj.ligand_mol_file)
        return out_data

    class Meta:
        model = models.SiteObservation
        fields = ("id", "vectors")


class MolpropsSerializer(serializers.ModelSerializer):
    props = serializers.SerializerMethodField()

    def get_props(self, obj):
        out_data = {}
        m = sanitize_mol(Chem.MolFromSmiles(obj.smiles))

        out_data["mol_log_p"] = Chem.Crippen.MolLogP(m)
        out_data["mol_wt"] = float(Chem.rdMolDescriptors.CalcExactMolWt(m))
        out_data["heavy_atom_count"] = Chem.Lipinski.HeavyAtomCount(m)
        out_data["heavy_atom_mol_wt"] = float(Descriptors.HeavyAtomMolWt(m))
        out_data["nhoh_count"] = Chem.Lipinski.NHOHCount(m)
        out_data["no_count"] = Chem.Lipinski.NOCount(m)
        out_data["num_h_acceptors"] = Chem.Lipinski.NumHAcceptors(m)
        out_data["num_h_donors"] = Chem.Lipinski.NumHDonors(m)
        out_data["num_het_atoms"] = Chem.Lipinski.NumHeteroatoms(m)
        out_data["num_rot_bonds"] = Chem.Lipinski.NumRotatableBonds(m)
        out_data["num_val_electrons"] = Descriptors.NumValenceElectrons(m)
        out_data["ring_count"] = Chem.Lipinski.RingCount(m)
        out_data["tpsa"] = Chem.rdMolDescriptors.CalcTPSA(m)

        return out_data

    class Meta:
        model = models.Compound
        fields = ("id", "props")


class GraphSerializer(serializers.ModelSerializer):
    graph = serializers.SerializerMethodField()
    graph_choice = settings.NEO4J_QUERY
    graph_auth = settings.NEO4J_AUTH

    def get_graph(self, obj):
        return get_full_graph(
            obj.smiles, self.graph_choice, self.graph_auth, isomericSmiles=True
        )

    class Meta:
        model = models.SiteObservation
        fields = ("id", "graph")


class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User
        fields = ('id', 'username', 'email', 'first_name', 'last_name')


# GET
class ActionTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ActionType
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
        sp_tags = models.SessionProjectTag.objects.filter(
            session_projects=obj.id
        ).values()
        return sp_tags

    class Meta:
        model = models.SessionProject
        fields = (
            'id',
            'title',
            'init_date',
            'description',
            'target',
            'project',
            'author',
            'tags',
            'session_project_tags',
        )


# (POST, PUT, PATCH)
class SessionProjectWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.SessionProject
        fields = '__all__'


# (GET, POST, PUT, PATCH)
class SessionActionsSerializer(serializers.ModelSerializer):
    actions = serializers.JSONField()

    class Meta:
        model = models.SessionActions
        fields = '__all__'


# GET
class SnapshotReadSerializer(serializers.ModelSerializer):
    author = UserSerializer()
    session_project = SessionProjectWriteSerializer()

    class Meta:
        model = models.Snapshot
        fields = (
            'id',
            'type',
            'title',
            'author',
            'description',
            'created',
            'data',
            'session_project',
            'parent',
            'children',
            'additional_info',
        )


# (POST, PUT, PATCH)
class SnapshotWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Snapshot
        fields = (
            'id',
            'type',
            'title',
            'author',
            'description',
            'created',
            'data',
            'session_project',
            'parent',
            'children',
            'additional_info',
        )


# (GET, POST, PUT, PATCH)
class SnapshotActionsSerializer(serializers.ModelSerializer):
    actions = serializers.JSONField()

    class Meta:
        model = models.SnapshotActions
        fields = '__all__'


class ComputedSetSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ComputedSet
        fields = '__all__'


class ComputedMoleculeSerializer(serializers.ModelSerializer):
    # performance issue
    # inspiration_frags = MoleculeSerializer(read_only=True, many=True)
    class Meta:
        model = models.ComputedMolecule
        fields = '__all__'


class ScoreDescriptionSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ScoreDescription
        fields = '__all__'


class NumericalScoreSerializer(serializers.ModelSerializer):
    score = ScoreDescriptionSerializer(read_only=True)

    class Meta:
        model = models.NumericalScoreValues
        fields = '__all__'


class TextScoreSerializer(serializers.ModelSerializer):
    score = ScoreDescriptionSerializer(read_only=True)

    class Meta:
        model = models.TextScoreValues
        fields = '__all__'


class ComputedMolAndScoreSerializer(serializers.ModelSerializer):
    numerical_scores = serializers.SerializerMethodField()
    text_scores = serializers.SerializerMethodField()
    pdb_info = serializers.SerializerMethodField()
    # score_descriptions = serializers.SerializerMethodField()

    class Meta:
        model = models.ComputedMolecule
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
        scores = models.NumericalScoreValues.objects.filter(compound=obj)
        score_dict = {}
        for score in scores:
            score_dict[score.score.name] = score.value
        return score_dict

    def get_text_scores(self, obj):
        scores = models.TextScoreValues.objects.filter(compound=obj)
        score_dict = {}
        for score in scores:
            score_dict[score.score.name] = score.value
        return score_dict

    def get_pdb_info(self, obj):
        # For this (new XCA) Fragalysis phase we do not support
        # PDB material in the ComputedMolecule. So instead of this (original code)
        # we now return a constant 'None'
        #        if obj.pdb:
        #            return obj.pdb.pdb_info.url
        #        else:
        #            return None

        # Unused arguments
        del obj

        return None

    # def get_score_descriptions(self, obj):
    #     descriptions = ScoreDescription.objects.filter(computed_set=obj.computed_set)
    #     desc_dict = {}
    #     for desc in descriptions:
    #         desc_dict[desc.name] = desc.description
    #     return desc_dict


# Class for customer Discourse API
class DiscoursePostWriteSerializer(serializers.Serializer):
    category_name = serializers.CharField(max_length=200)
    parent_category_name = serializers.CharField(
        max_length=200, initial=settings.DISCOURSE_PARENT_CATEGORY
    )
    category_colour = serializers.CharField(max_length=10, initial="0088CC")
    category_text_colour = serializers.CharField(max_length=10, initial="FFFFFF")
    post_title = serializers.CharField(max_length=200)
    post_content = serializers.CharField(max_length=2000)
    post_tags = serializers.JSONField()


class DictToCsvSerializer(serializers.Serializer):
    title = serializers.CharField(max_length=200)
    dict = serializers.DictField()


class TagCategorySerializer(serializers.ModelSerializer):
    class Meta:
        model = models.TagCategory
        fields = '__all__'


class SiteObservationTagSerializer(serializers.ModelSerializer):
    site_observations = serializers.PrimaryKeyRelatedField(
        many=True, queryset=models.SiteObservation.objects.all()
    )

    def create(self, validated_data):
        # populate 'upload_name' field at object creation
        validated_data['upload_name'] = validated_data['tag']
        return super().create(validated_data)

    class Meta:
        model = models.SiteObservationTag
        fields = '__all__'
        extra_kwargs = {
            "id": {"read_only": True},
            "upload_name": {"read_only": True},
            "tag_prefix": {"read_only": True},
        }


class SessionProjectTagSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.SessionProjectTag
        fields = '__all__'


# class FirstAssemblySerializer(serializers.ModelSerializer):
#     target = serializers.IntegerField()
#     pdb_file = serializers.SerializerMethodField()

#     def get_pdb_file(self, obj):
#         exp_upload = (
#             models.ExperimentUpload.objects.filter(
#                 target__id=obj.target,
#             )
#             .order_by('-commit_datetime')
#             .first()
#         )

#         yaml_path = (
#             Path(settings.MEDIA_ROOT)
#             .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
#             .joinpath(exp_upload.task_id)
#         )

#         # add unpacked zip directory
#         yaml_path = [d for d in list(yaml_path.glob("*")) if d.is_dir()][0]

#         # add upload_[d] dir
#         yaml_path = next(yaml_path.glob("upload_*"))

#         # last components of path, need for reconstruction later
#         comps = yaml_path.parts[-2:]

#         # and the file itself
#         yaml_path = yaml_path.joinpath(XTALFORMS_FILE)
#         logger.debug("assemblies path: %s", yaml_path)
#         if yaml_path.is_file():
#             with open(yaml_path, "r", encoding="utf-8") as file:
#                 contents = yaml.safe_load(file)
#                 try:
#                     assemblies = contents["assemblies"]
#                 except KeyError:
#                     logger.error("No 'assemblies' section in '%s'", XTALFORMS_FILE)
#                     return ''

#                 try:
#                     first = list(assemblies.values())[0]
#                 except IndexError:
#                     logger.error("No assemblies in 'assemblies' section")
#                     return ''

#                 try:
#                     reference = first["reference"]
#                 except KeyError:
#                     logger.error("No assemblies in 'assemblies' section")
#                     return ''

#                 ref_path = (
#                     Path(settings.TARGET_LOADER_MEDIA_DIRECTORY)
#                     .joinpath(exp_upload.task_id)
#                     .joinpath(comps[0])
#                     .joinpath(comps[1])
#                     .joinpath("crystallographic_files")
#                     .joinpath(reference)
#                     .joinpath(f"{reference}.pdb")
#                 )
#                 logger.debug('ref_path: %s', ref_path)
#                 if Path(settings.MEDIA_ROOT).joinpath(ref_path).is_file():
#                     request = self.context.get('request', None)
#                     if request is not None:
#                         return request.build_absolute_uri(
#                             Path(settings.MEDIA_URL).joinpath(ref_path)
#                         )
#                     else:
#                         return ''
#                 else:
#                     logger.error("Reference pdb file doesn't exist")
#                     return ''
#         else:
#             logger.error("'%s' missing", XTALFORMS_FILE)
#             return ''

#     class Meta:
#         model = models.QuatAssembly
#         fields = ('target', 'pdb_file')


# class TargetMoleculesSerializer(serializers.ModelSerializer):
#     template_protein = serializers.SerializerMethodField()
#     zip_archive = serializers.SerializerMethodField()
#     metadata = serializers.SerializerMethodField()
#     sequences = serializers.SerializerMethodField()
#     molecules = serializers.SerializerMethodField()
#     tags_info = serializers.SerializerMethodField()
#     tag_categories = serializers.SerializerMethodField()

#     def get_template_protein(self, obj):
#         return template_protein(obj)

#     def get_zip_archive(self, obj):
#         # The if-check is because the filefield in target has null=True.
#         # Note that this link will not work on local
#         if hasattr(obj, 'zip_archive') and obj.zip_archive.name:
#             return urljoin(get_https_host(self.context["request"]), obj.zip_archive.url)
#         return

#     def get_metadata(self, obj):
#         if hasattr(obj, 'metadata') and obj.metadata.name:
#             return urljoin(get_https_host(self.context["request"]), obj.metadata.url)
#         return

#     def get_sequences(self, obj):
#         return protein_sequences(obj)

#     def get_molecules(self, obj):
#         mols = models.SiteObservation.objects.filter(
#             experiment__experiment_upload__target__id=obj.id
#         ).annotate(
#             # NB! some of the fields are just renamed here, avoiding
#             # that would simplify things here and remove some lines of
#             # code. but that means front-end code needs to know about
#             # the changes
#             protein_code=F('code'),
#             molecule_protein=F('experiment__pdb_info'),
#             sdf_info=F('ligand_mol_file'),
#             lig_id=F('seq_id'),
#             tags_set=ArrayAgg("siteobservationtag__pk"),
#         )
#         fields = [
#             "id",
#             "smiles",
#             "cmpd",
#             "protein_code",
#             "molecule_protein",
#             "lig_id",
#             "chain_id",
#             "sdf_info",
#             "tags_set",
#         ]

#         logger.debug("%s", mols)
#         logger.debug("%s", mols.values(*fields))

#         molecules = [
#             {'data': k, 'tags_set': k['tags_set']} for k in mols.values(*fields)
#         ]

#         return molecules

#     def get_tags_info(self, obj):
#         tags = models.SiteObservationTag.objects.filter(target_id=obj.id)
#         tags_info = []
#         for tag in tags:
#             tag_data = models.SiteObservationTag.objects.filter(id=tag.id).values()
#             tag_coords = SiteObservationGroup.objects.filter(
#                 id=tag.mol_group_id
#             ).values('x_com', 'y_com', 'z_com')
#             tag_dict = {'data': tag_data, 'coords': tag_coords}
#             tags_info.append(tag_dict)

#         return tags_info

#     def get_tag_categories(self, obj):
#         tag_categories = (
#             models.TagCategory.objects.filter(siteobservationtag__target_id=obj.id)
#             .distinct()
#             .values()
#         )
#         return tag_categories

#     class Meta:
#         model = models.Target
#         fields = (
#             "id",
#             "title",
#             "project_id",
#             "default_squonk_project",
#             "template_protein",
#             "metadata",
#             "zip_archive",
#             "upload_status",
#             "sequences",
#             "molecules",
#             "tags_info",
#             "tag_categories",
#         )


class DownloadStructuresSerializer(serializers.Serializer):
    target_name = serializers.CharField(max_length=200, default=None, allow_blank=True)
    proteins = serializers.CharField(max_length=5000, default='', allow_blank=True)
    all_aligned_structures = serializers.BooleanField(default=False)
    pdb_info = serializers.BooleanField(default=False)
    cif_info = serializers.BooleanField(default=False)
    mtz_info = serializers.BooleanField(default=False)
    diff_file = serializers.BooleanField(default=False)
    event_file = serializers.BooleanField(default=False)
    sigmaa_file = serializers.BooleanField(default=False)
    map_info = serializers.BooleanField(default=False)
    single_sdf_file = serializers.BooleanField(default=False)
    metadata_info = serializers.BooleanField(default=False)
    static_link = serializers.BooleanField(default=False)
    file_url = serializers.CharField(max_length=200, default='', allow_blank=True)
    trans_matrix_info = serializers.BooleanField(default=False)


# Start of Serializers for Squonk Jobs
# (GET)
class JobFileTransferReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobFileTransfer
        fields = '__all__'


# (POST, PUT, PATCH)
class JobFileTransferWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobFileTransfer
        fields = ("snapshot", "target", "squonk_project", "proteins", "compounds")


# (GET)
class JobRequestReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobRequest
        fields = '__all__'


# (POST, PUT, PATCH)
class JobRequestWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobRequest
        fields = (
            "squonk_job_name",
            "snapshot",
            "target",
            "squonk_project",
            "squonk_job_spec",
        )


class JobCallBackReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobRequest
        fields = (
            "job_status",
            "job_status_datetime",
            "squonk_job_name",
            "squonk_job_spec",
            "upload_status",
        )


class JobCallBackWriteSerializer(serializers.ModelSerializer):
    SQUONK_STATUS = ['PENDING', 'STARTED', 'SUCCESS', 'FAILURE', 'RETRY', 'REVOKED']
    job_status = serializers.ChoiceField(choices=SQUONK_STATUS, default="PENDING")
    state_transition_time = serializers.DateTimeField(source='job_status_datetime')

    class Meta:
        model = models.JobRequest
        fields = ("job_status", "state_transition_time")


class TargetExperimentReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.ExperimentUpload
        fields = '__all__'


class TargetExperimentWriteSerializer(serializers.ModelSerializer):
    target_access_string = serializers.CharField(label='Target Access String')
    contact_email = serializers.EmailField(required=False, default=None)

    def validate(self, data):
        """Verify TAS is correctly formed."""
        success, error_msg = validate_tas(data['target_access_string'])
        if not success:
            raise serializers.ValidationError({"target_access_string": error_msg})
        return data

    class Meta:
        model = models.ExperimentUpload
        fields = (
            'target_access_string',
            'contact_email',
            'file',
        )


class TargetExperimentDownloadSerializer(serializers.ModelSerializer):
    filename = serializers.CharField()

    class Meta:
        model = models.ExperimentUpload
        fields = (
            'project',
            'target',
            'filename',
        )


class JobOverrideReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobOverride
        fields = '__all__'


class JobOverrideWriteSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.JobOverride
        fields = ('override',)


class SiteObservationReadSerializer(serializers.ModelSerializer):
    compound_code = serializers.StringRelatedField()
    prefix_tooltip = serializers.StringRelatedField()

    class Meta:
        model = models.SiteObservation
        fields = '__all__'


class CanonSiteReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CanonSite
        fields = '__all__'


class CanonSiteConfReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CanonSiteConf
        fields = '__all__'


class XtalformSiteReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.XtalformSite
        fields = '__all__'
