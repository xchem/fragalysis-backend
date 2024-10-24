import contextlib
import logging
from pathlib import Path
from urllib.parse import urljoin

import yaml
from django.conf import settings
from django.contrib.auth.models import User
from django.db import IntegrityError, transaction
from django.db.models import Count
from frag.network.decorate import get_3d_vects_for_mol, get_vect_indices_for_mol
from frag.network.query import get_full_graph
from rdkit import Chem
from rdkit.Chem import Descriptors
from rest_framework import serializers
from rest_framework.exceptions import PermissionDenied

from api.security import ISPyBSafeQuerySet
from api.utils import draw_mol, validate_tas
from viewer import models
from viewer.target_loader import XTALFORMS_FILE
from viewer.target_set_upload import sanitize_mol
from viewer.utils import get_https_host

logger = logging.getLogger(__name__)

_ISPYB_SAFE_QUERY_SET = ISPyBSafeQuerySet()


class ValidateProjectMixin:
    """Mixin for serializers to check if user is allowed to create objects.

    Requires a 'filter_permissions' member in the corresponding View.
    This is used to navigate to the Project object from the data map
    given to the validate() method.
    """

    def validate(self, data):
        # User must be logged in
        user = self.context['request'].user  # type: ignore [attr-defined]
        if not user or not user.is_authenticated:
            raise serializers.ValidationError("You must be logged in")
        view = self.context['view']  # type: ignore [attr-defined]
        if not hasattr(view, "filter_permissions"):
            raise AttributeError(
                "The view object must define a 'filter_permissions' property"
            )

        # We expect a filter_permissions string (defined in the View) like this...
        #     "compound__project_id"
        # In this example the supplied data map is therefore expected to have a
        # "compound" key (which we extract into a variable called 'base_object_key').
        # We use the 2nd half of the string (which we call 'project_path')
        # to get to the Project object from 'data["compound"]'.
        #
        # If the filter_permissions string has no 2nd half (e.g. it's simply 'project_id')
        # then the data is clearly expected to contain the Project object itself.

        base_object_key, project_path = view.filter_permissions.split('__', 1)
        base_start_obj = data[base_object_key]
        # Assume we're using the base object,
        # but swap it out of there's a project path.
        project_obj = base_start_obj
        if project_path:
            try:
                project_obj = getattr(base_start_obj, project_path)
            except AttributeError as exc:
                # Something's gone wrong trying to lookup the project.
                # Log some 'interesting' contextual information...
                logger.info('context=%s', self.context)  # type: ignore [attr-defined]
                logger.info('data=%s', data)
                logger.info('view=%s', view.__class__.__name__)
                logger.info('view.filter_permissions=%s', view.filter_permissions)
                # Get the object's content and dump it for analysis...
                bso_class_name = base_start_obj.__class__.__name__
                msg = f"There is no Project at '{project_path}' ({view.filter_permissions})"
                logger.error(
                    "%s - base_start_obj=%s vars(base_start_obj)=%s",
                    msg,
                    bso_class_name,
                    vars(base_start_obj),
                )
                raise serializers.ValidationError(msg) from exc
        assert project_obj
        # Now get the proposals from the Project(s)...
        object_proposals = project_obj.values_list('title', flat=True)
        if not object_proposals:
            raise PermissionDenied(
                detail="Authority cannot be granted - the object is not a part of any Project"
            )

        # Now we have the proposals (Project titles) the object belongs to,
        # has the user been associated (in IPSpyB) with any of them?
        # We can always see (GET) objects that are open to the public.
        restrict_public = False if self.context['request'].method == 'GET' else True  # type: ignore [attr-defined]
        if not _ISPYB_SAFE_QUERY_SET.user_is_member_of_any_given_proposals(
            user=user,
            proposals=object_proposals,
            restrict_public_to_membership=restrict_public,
        ):
            raise PermissionDenied(
                detail="Your authority to access this object has not been given"
            )

        # OK if we get here...
        return data


class FileSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.File
        fields = "__all__"


class CompoundIdentifierTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CompoundIdentifierType
        fields = '__all__'


class CompoundIdentifierSerializer(ValidateProjectMixin, serializers.ModelSerializer):
    class Meta:
        model = models.CompoundIdentifier
        fields = '__all__'


class TargetSerializer(serializers.ModelSerializer):
    template_protein = serializers.SerializerMethodField()
    zip_archive = serializers.SerializerMethodField()
    metadata = serializers.SerializerMethodField()

    def get_template_protein_path(self, experiment_upload) -> Path | None:
        yaml_path = experiment_upload.get_upload_path()

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
                    return None

                try:
                    first = list(assemblies.values())[0]
                except IndexError:
                    logger.error("No assemblies in 'assemblies' section")
                    return None

                try:
                    reference = first["reference"]
                except KeyError:
                    logger.error("No assemblies in 'assemblies' section")
                    return None

                ref_path = (
                    Path(settings.TARGET_LOADER_MEDIA_DIRECTORY)
                    .joinpath(experiment_upload.target.zip_archive.name)
                    .joinpath(experiment_upload.upload_data_dir)
                    .joinpath("crystallographic_files")
                    .joinpath(reference)
                    .joinpath(f"{reference}.pdb")
                )
                logger.debug('ref_path: %s', ref_path)
                if Path(settings.MEDIA_ROOT).joinpath(ref_path).is_file():
                    return ref_path
                else:
                    logger.error("Reference pdb file doesn't exist")
                    return None
        else:
            logger.error("'%s' missing", XTALFORMS_FILE)
            return None

    def get_template_protein(self, obj):
        # loop through exp uploads from latest to earliest, and try to
        # find template protein
        for exp_upload in models.ExperimentUpload.objects.filter(
            target=obj,
        ).order_by('-commit_datetime'):
            path = self.get_template_protein_path(exp_upload)
            if path is None:
                continue
            else:
                request = self.context.get('request', None)
                if request is not None:
                    return request.build_absolute_uri(
                        Path(settings.MEDIA_URL).joinpath(path)
                    )
                else:
                    return None
        return None

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

    class Meta:
        model = models.Target
        fields = (
            "id",
            "title",
            "display_name",
            "project",
            "default_squonk_project",
            "template_protein",
            "metadata",
            "zip_archive",
            "upload_status",
        )
        extra_kwargs = {
            "id": {"read_only": True},
            "title": {"read_only": True},
            "project": {"read_only": True},
            "default_squonk_project": {"read_only": True},
            "template_protein": {"read_only": True},
            "metadata": {"read_only": True},
            "zip_archive": {"read_only": True},
            "upload_status": {"read_only": True},
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


class ComputedSetSerializer(ValidateProjectMixin, serializers.ModelSerializer):
    class Meta:
        model = models.ComputedSet
        fields = '__all__'


class ComputedSetDownloadSerializer(serializers.ModelSerializer):
    # validation is not called, so no reason to use it
    class Meta:
        model = models.ComputedSet
        fields = ('name',)


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
    filename = serializers.CharField()
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


class DownloadStructuresSerializer(serializers.Serializer):
    target_name = serializers.CharField(max_length=200, default=None, allow_blank=True)
    target_access_string = serializers.CharField(
        max_length=200, default=None, allow_blank=True
    )
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
class JobFileTransferReadSerializer(ValidateProjectMixin, serializers.ModelSerializer):
    class Meta:
        model = models.JobFileTransfer
        fields = '__all__'


# (POST, PUT, PATCH)
class JobFileTransferWriteSerializer(ValidateProjectMixin, serializers.ModelSerializer):
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


class TargetExperimentReadSerializer(ValidateProjectMixin, serializers.ModelSerializer):
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

    ligand_mol_file = serializers.SerializerMethodField()

    def get_ligand_mol_file(self, obj):
        contents = ''
        if obj.ligand_mol:
            path = Path(settings.MEDIA_ROOT).joinpath(obj.ligand_mol.name)
            with contextlib.suppress(TypeError, FileNotFoundError):
                with open(path, "r", encoding="utf-8") as f:
                    contents = f.read()

        return contents

    class Meta:
        model = models.SiteObservation
        fields = '__all__'


class CanonSiteReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CanonSite
        fields = '__all__'


class ExperimentReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Experiment
        fields = '__all__'


class CanonSiteConfReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.CanonSiteConf
        fields = '__all__'


class XtalformSiteReadSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.XtalformSite
        fields = '__all__'


class PoseSerializer(ValidateProjectMixin, serializers.ModelSerializer):
    site_observations = serializers.PrimaryKeyRelatedField(
        many=True,
        queryset=models.SiteObservation.objects.all(),
    )
    main_site_observation = serializers.PrimaryKeyRelatedField(
        required=False, default=None, queryset=models.SiteObservation.objects.all()
    )
    main_site_observation_cmpd_code = serializers.SerializerMethodField()

    def get_main_site_observation_cmpd_code(self, obj):
        return obj.main_site_observation.cmpd.compound_code

    def update_pose_main_observation(self, main_obvs):
        if main_obvs.pose.site_observations.count() == 1:
            main_obvs.pose.delete()

    def remove_empty_poses(self, pose):
        # fmt: off
        empty_poses = models.Pose.objects.filter(
            canon_site=pose.canon_site,
            compound=pose.compound,
        ).annotate(
            num_obvs=Count('site_observations'),
        ).filter(
            num_obvs=0,
        )
        # fmt: on
        logger.debug('empty_poses: %s', empty_poses)
        empty_poses.delete()

    def create(self, validated_data):
        logger.debug('validated_data: %s', validated_data)
        try:
            with transaction.atomic():
                # only situation with suggested main being a main in
                # other pose is when it's the only observation in
                # pose, the other cases must have been caught by
                # validation
                self.update_pose_main_observation(
                    validated_data['main_site_observation']
                )
                pose = super().create(validated_data)
                self.remove_empty_poses(pose)
                return pose

        except IntegrityError as exc:
            raise serializers.ValidationError({'errors': exc.args[0]})

    def update(self, instance, validated_data):
        logger.debug('validated_data: %s', validated_data)

        try:
            with transaction.atomic():
                # TODO: func here to release the main
                pose = super().update(instance, validated_data)
                pose.save()
                self.remove_empty_poses(pose)
                return pose

        except IntegrityError as exc:
            raise serializers.ValidationError({'errors': exc.args[0]})

    def validate(self, data):
        """Validate data sent to create or update pose.

        First look at the main_site_observation. If this is not set,
        try to find a suitable one from the list of observations. If
        is given, check if it already isn't a main observation to some
        othe pose.

        Then check observation and make sure that:
        - they have the same canon_site and compound as the pose being created
        - the pose they're being removed from is deleted when empty

        """
        logger.info('+ validate data: %s', data)

        data = super().validate(data)

        template = (
            "Site observation {} cannot be assigned to pose because "
            + "pose's {} ({}) doesn't match with observation's ({})"
        )

        messages: dict[str, list[str]] = {
            'site_observations': [],
            'main_site_observation': [],
        }
        validated_observations = []

        # populate data struct with missing attributes but only when
        # they're missing, not simply undefined. not sure why the
        # serializer isn't doing this here
        instance = getattr(self, 'instance', None)
        if 'site_observations' not in data.keys():
            data['site_observations'] = (
                instance.site_observations.all() if instance else []
            )
        if 'main_site_observation' not in data.keys():
            data['main_site_observation'] = (
                instance.main_site_observation if instance else None
            )
        if 'compound' not in data.keys():
            data['compound'] = instance.compound if instance else None
        if 'canon_site' not in data.keys():
            data['canon_site'] = instance.canon_site if instance else None

        if not data['site_observations']:
            if hasattr(self, 'instance'):
                messages['site_observations'].append(
                    'Cannot remove all site observations from pose'
                )
            else:
                messages['site_observations'].append('Cannot create empty pose')

        if data['main_site_observation']:
            logger.debug('main observation given, trying to set one')
            main_pose = getattr(data['main_site_observation'], 'main_pose', None)
            if (
                main_pose
                and main_pose != instance
                and main_pose.site_observations.count() > 1
            ):
                # checking length, because if single observation in
                # pose, it would leave an empty pose which will be
                # deleted. otherwise don't allow setting main
                messages['main_site_observation'].append(
                    'Requested main_site_observation already the main observation '
                    + f'in pose {main_pose.display_name}'
                )

        for so in data['site_observations']:
            all_good = True

            # only try to set one if not given in request data
            if not data['main_site_observation']:
                logger.debug('no main observation given, trying to find one')
                if not hasattr(so, 'main_pose'):
                    data['main_site_observation'] = so

            logger.debug('processing observation: %s', so)
            if so.canon_site_conf.canon_site != data['canon_site']:
                all_good = False
                messages['site_observations'].append(
                    template.format(
                        so.code,
                        'canonical site',
                        data['canon_site'].name,
                        so.canon_site_conf.canon_site.name,
                    )
                )
            if so.cmpd != data['compound']:
                all_good = False
                messages['site_observations'].append(
                    template.format(
                        so.code,
                        'compound',
                        data['compound'].compound_code,
                        so.cmpd.compound_code,
                    )
                )
            logger.debug('accumulated messages: %s', messages)
            if all_good:
                validated_observations.append(so)

        # did we get main observaton?
        if not data['main_site_observation']:
            messages['main_site_observation'].append(
                'No suitable observation found to set as main_site_observation'
            )

        if any(messages.values()):
            raise serializers.ValidationError(messages)
        else:
            data['site_observations'] = validated_observations
            return data

    class Meta:
        model = models.Pose
        fields = (
            'id',
            'display_name',
            'canon_site',
            'compound',
            'main_site_observation',
            'site_observations',
            'main_site_observation_cmpd_code',
        )


class MetadataUploadSerializer(serializers.Serializer):
    filename = serializers.FileField()
    target = serializers.CharField()
    target_access_string = serializers.CharField()
