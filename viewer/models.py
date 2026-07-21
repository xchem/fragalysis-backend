import contextlib
import logging
import os
import uuid
from dataclasses import dataclass
from pathlib import Path

from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.core.serializers.json import DjangoJSONEncoder
from django.core.validators import MinLengthValidator
from django.db import IntegrityError, models, transaction
from django.utils import timezone
from pgvector.django import HalfVectorField, HnswIndex
from shortuuid.django_fields import ShortUUIDField
from simple_history.models import HistoricalRecords

from .managers import (
    AssayResultDataManager,
    CanonSiteConfDataManager,
    CanonSiteDataManager,
    CompoundDataManager,
    CompoundIdentifierDataManager,
    ComputedSetDataManager,
    ExperimentDataManager,
    ExperimentUploadDataManager,
    PoseDataManager,
    QuatAssemblyDataManager,
    ResultUploadDataManager,
    SessionActionsDataManager,
    SiteObservationDataManager,
    SiteObservationQualityStatusDataManager,
    SnapshotActionsDataManager,
    SnapshotDataManager,
    XtalformDataManager,
    XtalformQuatAssemblyDataManager,
    XtalformSiteDataManager,
)

logger = logging.getLogger(__name__)


@dataclass
class Vector3d:
    start_x: float
    start_y: float
    start_z: float
    end_x: float
    end_y: float
    end_z: float
    number: int
    vector_type: str
    smiles: str
    site_observation: int
    cmpd_id: int


class Project(models.Model):
    title = models.TextField(null=False, unique=True)
    alias = models.TextField(null=True)
    init_date = models.DateTimeField(auto_now_add=True)
    user_id = models.ManyToManyField(User)
    open_to_public = models.BooleanField(default=False)

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<Project %r %r %r>" % (self.id, self.title, self.open_to_public)


class UserRole(models.Model):
    LOADER_ROLE = "Loader"

    name = models.TextField(unique=True)
    users = models.ManyToManyField(User, related_name="roles", blank=True)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<UserRole %r %r>" % (self.id, self.name)


class Target(models.Model):
    PENDING = "PENDING"
    STARTED = "STARTED"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"
    RETRY = "RETRY"
    REVOKED = "REVOKED"
    STATUS = (
        (PENDING, 'PENDING'),  # Initial state when queued
        (STARTED, 'STARTED'),  # File transfer started
        (SUCCESS, 'SUCCESS'),  # File transfer finished successfully
        (FAILURE, 'FAILURE'),  # File transfer failed
        (RETRY, 'RETRY'),
        (REVOKED, 'REVOKED'),
    )

    title = models.CharField(max_length=200, help_text="A title, i.e. Mpro")
    display_name = models.TextField(null=False, blank=True)

    init_date = models.DateTimeField(auto_now_add=True)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)

    uniprot_id = models.CharField(
        max_length=100,
        null=True,
        help_text="The uniprot ID id for the target. A unique key",
    )
    metadata = models.FileField(
        upload_to="metadata/",
        null=True,
        max_length=255,
        help_text='Optional file upload defining metadata about the target.'
        ' Can be used to add custom site labels',
    )
    zip_archive = models.FileField(
        upload_to="archive/",
        null=True,
        max_length=255,
        help_text='Link to zip file created from targets uploaded with the loader',
    )
    default_squonk_project = models.CharField(max_length=200, null=True)
    upload_task_id = models.CharField(
        null=True, max_length=50, help_text='The Task ID of upload Celery task)'
    )
    upload_status = models.CharField(
        choices=STATUS,
        null=True,
        max_length=7,
        help_text='Identifies the status of the upload.'
        ' Will only be updated at the end of the process',
    )
    upload_progress = models.DecimalField(
        null=True,
        max_digits=5,
        decimal_places=2,
        help_text='Intended to be used as an indication of progress (0 to 100%)',
    )
    upload_datetime = models.DateTimeField(
        null=True, help_text='The datetime the upload was completed'
    )
    # this is to be deprecated and used as part of the settings
    alias_order = ArrayField(models.TextField(), null=True)
    short_name = models.TextField(null=True, blank=True)
    long_name = models.TextField(null=True, blank=True)
    organism = models.TextField(null=True, blank=True)
    external_url = models.URLField(null=True, blank=True)
    external_url_display_name = models.TextField(null=True, blank=True)
    settings = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        blank=True,
    )

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<Target %r %r %r %r>" % (
            self.id,
            self.title,
            self.display_name,
            self.project,
        )

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "title",
                    "project",
                ],
                name="unique_target_in_project",
            ),
        ]


class ExperimentUpload(models.Model):
    LOADING = "LOADING"
    LOADED = "LOADED"
    FAILURE = "FAILURE"
    STATUS = (
        (LOADING, "LOADING"),
        (LOADED, "LOADED"),
        (FAILURE, "FAILURE"),
    )
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    file = models.FileField(upload_to="experiment-upload/", max_length=255)
    commit_datetime = models.DateTimeField(
        help_text="The UTC datetime the upload was committed"
    )
    committer = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        help_text="The user committing the original file."
        " This user may not be the author of the file",
    )
    complete_datetime = models.DateTimeField(
        null=True,
        help_text="The UTC datetime the upload finished."
        " It can be considered a success"
        " if the status is LOADED",
    )
    task_id = models.CharField(
        null=True, max_length=50, help_text="Celery task ID responsible for the upload"
    )
    status = models.CharField(choices=STATUS, default=LOADING, max_length=7)
    message = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        blank=True,
        help_text="Any message or task info associated with the upload."
        " Used for upload audit trail",
    )
    neighbourhood_transforms = models.FileField(
        upload_to="experiment-upload/", max_length=255, null=True
    )
    conformer_site_transforms = models.FileField(
        upload_to="experiment-upload/",
        max_length=255,
        null=True,
    )
    reference_structure_transforms = models.FileField(
        upload_to="experiment-upload/",
        max_length=255,
        null=True,
    )
    assembly_transforms = models.FileField(
        upload_to="experiment-upload/",
        max_length=255,
        null=True,
    )
    upload_data_dir = models.TextField(null=True)
    upload_version = models.PositiveSmallIntegerField(default=1)
    data_version_major = models.PositiveSmallIntegerField(default=0)
    data_version_minor = models.PositiveSmallIntegerField(default=0)

    objects = models.Manager()
    upload_manager = ExperimentUploadDataManager()

    def __str__(self) -> str:
        return f"{self.project}"

    def __repr__(self) -> str:
        return "<ExperimentUpload %r %r %r>" % (self.id, self.project, self.target)

    def get_upload_path(self):
        return (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(self.target.zip_archive.name)
            .joinpath(self.upload_data_dir)
        )

    def get_download_path(self):
        """The path to the original uploaded file, used during downloads"""
        return (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(Path(str(self.file)))
        )


class QualityStatusType(models.Model):
    status = models.TextField(primary_key=True)


class RefinementStatusType(models.Model):
    code = models.IntegerField(blank=True)
    description = models.TextField(blank=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "code",
                ],
                name="unique_refinement_code",
            ),
        ]

    def __repr__(self) -> str:
        return "<RefinementStatusType %r %r %r>" % (
            self.id,
            self.code,
            self.description,
        )


class ExperimentStatusType(models.Model):
    status_code = models.IntegerField(blank=True, null=True)
    status = models.TextField(null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "status_code",
                ],
                name="unique_experiment_status",
            ),
        ]


class Experiment(models.Model):
    experiment_upload = models.ForeignKey(ExperimentUpload, on_delete=models.CASCADE)
    code = models.TextField(null=True)
    status = models.ForeignKey(
        ExperimentStatusType,
        to_field='status_code',
        on_delete=models.SET_NULL,
        null=True,
    )
    pdb_info = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    pdb_info_source_file = models.TextField(null=True)
    mtz_info = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    mtz_info_source_file = models.TextField(null=True)
    cif_info = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    cif_info_source_file = models.TextField(null=True)
    map_info = ArrayField(models.FileField(max_length=255), null=True)
    map_info_source_files = ArrayField(models.TextField(null=True), null=True)
    type = models.PositiveSmallIntegerField(null=True)
    pdb_sha256 = models.TextField(null=True)
    prefix_tooltip = models.TextField(null=True)
    code_prefix = models.TextField(null=True)
    compounds = models.ManyToManyField(
        "Compound",
        through="ExperimentCompound",
        through_fields=("experiment", "compound"),
    )
    # need to set null=True due to the data saving order
    xtalform = models.ForeignKey("Xtalform", null=True, on_delete=models.CASCADE)
    refinement_outcome = models.ForeignKey(
        RefinementStatusType,
        on_delete=models.SET_NULL,
        null=True,
    )

    # new fields from SoakDb (issue 1147)
    cchalf_high_res_shell = models.FloatField(null=True, blank=True)
    cchalf_overall = models.FloatField(null=True, blank=True)
    completeness_high_res_shell = models.FloatField(null=True, blank=True)
    completeness_overall = models.FloatField(null=True, blank=True)
    # potential for lookup? seems very standardised
    crystal_mounting_result = models.TextField(null=True, blank=True)
    data_collection_date = models.DateTimeField(null=True)
    data_collection_outcome = models.TextField(null=True, blank=True)
    # same as code
    # dataset = models.TextField(null=True, blank=True)
    date_model_last_updated = models.DateTimeField(null=True)
    date_status_updated = models.DateTimeField(null=True)
    date_refined = models.DateTimeField(null=True)
    dimple_rfree = models.FloatField(null=True, blank=True)
    dimple_rwork = models.FloatField(null=True, blank=True)
    experiment_comments = models.TextField(null=True, blank=True)
    experiment_status = models.TextField(null=True, blank=True)
    experiment_type = models.TextField(null=True, blank=True)
    experiment_start_date = models.DateTimeField(null=True)
    final_compound_concentration_mm = models.FloatField(null=True, blank=True)
    high_resolution = models.FloatField(null=True, blank=True)
    isig_i_overall = models.FloatField(null=True, blank=True)
    isig_i_high_res_shell = models.FloatField(null=True, blank=True)
    library = models.TextField(null=True, blank=True)
    library_plate = models.TextField(null=True, blank=True)
    ligand_confidence = models.TextField(null=True, blank=True)
    ligand_correlation_coefficient = models.TextField(null=True, blank=True)
    model_last_updated_by = models.ForeignKey(
        User,
        null=True,
        on_delete=models.CASCADE,
        related_name="+",
    )
    modelled_smiles = models.TextField(null=True, blank=True)
    panddarun = models.TextField(null=True, blank=True)
    pdb_code = models.TextField(null=True, blank=True)
    processing_pipeline = models.TextField(null=True, blank=True)
    refined_by = models.ForeignKey(
        User,
        null=True,
        on_delete=models.CASCADE,
        related_name="+",
    )
    refinement_comment = models.TextField(null=True, blank=True)
    refinement_rfree = models.FloatField(null=True, blank=True)
    refinement_rwork = models.FloatField(null=True, blank=True)
    soakdb_entry = models.TextField(null=True, blank=True)
    soaking_time = models.DurationField(null=True, blank=True)
    source_well = models.TextField(null=True, blank=True)
    space_group = models.TextField(null=True, blank=True)
    unit_cell_dimensions = ArrayField(models.FloatField(), null=True)
    refinement_resolution = models.FloatField(null=True, blank=True)

    objects = models.Manager()
    filter_manager = ExperimentDataManager()

    def __str__(self) -> str:
        return f"{self.code}"

    def __repr__(self) -> str:
        return "<Experiment %r %r %r>" % (self.id, self.code, self.experiment_upload)


class Compound(models.Model):
    """Information about a compound, which is a unique 2D molecule"""

    inchi = models.TextField(unique=False, db_index=True)
    smiles = models.CharField(max_length=255, db_index=True)
    # rdkit representation of smiles field for structure-based
    # search. Internally rdkit mol type
    smiles_mol = models.TextField(editable=False, null=True)
    compound_code = models.TextField(null=True)
    current_identifier = models.OneToOneField(
        'CompoundIdentifier',
        blank=True,
        null=True,
        on_delete=models.SET_NULL,
        related_name='+',
        help_text='The preferred alias for this compound.',
    )
    project_id = models.ManyToManyField(Project)
    inspirations = models.ManyToManyField(
        "SiteObservation",
        blank=True,
        help_text='Foreign key link to any number of 3D Molecules that inspired'
        ' the design of this compound',
    )
    description = models.TextField(blank=True, null=True)
    comments = models.TextField(blank=True, null=True)
    inchi_key = models.CharField(db_index=True, max_length=27, blank=True)
    ligand_name = models.TextField(blank=True, default='LIG')
    modeled_smiles_soakdb = models.TextField(blank=True, null=True)
    modeled_smiles_canon = models.TextField(blank=True, null=True)
    soaked_smiles_soakdb = models.TextField(blank=True, null=True)
    soaked_smiles_canon = models.TextField(blank=True, null=True)

    objects = models.Manager()
    filter_manager = CompoundDataManager()

    def __str__(self) -> str:
        return f"{self.smiles}"

    def __repr__(self) -> str:
        return "<Compound %r %r %r>" % (self.id, self.smiles, self.inchi)


class ExperimentCompound(models.Model):
    experiment = models.ForeignKey(
        Experiment,
        null=False,
        on_delete=models.CASCADE,
    )
    compound = models.ForeignKey(
        Compound,
        null=False,
        on_delete=models.CASCADE,
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "experiment",
                    "compound",
                ],
                name="unique_experimentcompound",
            ),
        ]


class QuatAssembly(models.Model):
    chains = models.TextField()
    name = models.TextField()
    assembly_num = models.IntegerField(
        null=True, help_text="numeric assembly id (enumerated on creation)"
    )

    objects = models.Manager()
    filter_manager = QuatAssemblyDataManager()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<QuatAssembly %r %r %r>" % (self.id, self.name, self.chains)


class Xtalform(models.Model):
    name = models.TextField(null=True)
    space_group = models.TextField(null=True)
    unit_cell_info = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    quat_assembly = models.ManyToManyField(
        QuatAssembly,
        through="XtalformQuatAssembly",
        through_fields=("xtalform", "quat_assembly"),
    )
    xtalform_num = models.IntegerField(
        null=True, help_text="numeric xtalform id (enumerated on creation)"
    )

    objects = models.Manager()
    filter_manager = XtalformDataManager()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Xtalform %r %r>" % (self.id, self.name)


class XtalformQuatAssembly(models.Model):
    xtalform = models.ForeignKey(
        Xtalform,
        null=False,
        on_delete=models.CASCADE,
    )
    quat_assembly = models.ForeignKey(
        QuatAssembly,
        null=False,
        on_delete=models.CASCADE,
    )
    assembly_id = models.TextField()
    chains = models.TextField()

    objects = models.Manager()
    filter_manager = XtalformQuatAssemblyDataManager()

    def __str__(self) -> str:
        return f"XtalformQuatAssembly {self.xtalform} {self.quat_assembly} {self.assembly_id}"

    def __repr__(self) -> str:
        return "<XtalformQuatAssembly %r %r %r %r>" % (
            self.id,
            self.xtalform,
            self.quat_assembly,
            self.assembly_id,
        )

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "xtalform",
                    "quat_assembly",
                    "assembly_id",
                ],
                name="unique_xtalformquatassembly",
            ),
        ]


class Versionable(models.Model):
    superseded = models.BooleanField(null=False, default=False)
    version = models.PositiveSmallIntegerField(null=False, default=1)

    class Meta:
        abstract = True


class CanonSite(Versionable, models.Model):
    name = models.TextField()
    residues = models.JSONField(encoder=DjangoJSONEncoder)
    # TODO: missing in db, check if correct, (might be correct, but might not)
    ref_conf_site = models.OneToOneField(
        "CanonSiteConf", null=True, on_delete=models.CASCADE
    )
    canon_site_num = models.IntegerField(
        null=True, help_text="numeric canon site id (enumerated on creation)"
    )
    centroid_res = models.TextField(null=True)

    objects = models.Manager()
    filter_manager = CanonSiteDataManager()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CanonSite %r %r>" % (self.id, self.name)


class XtalformSite(Versionable, models.Model):
    xtalform = models.ForeignKey(Xtalform, on_delete=models.CASCADE)
    canon_site = models.ForeignKey(CanonSite, on_delete=models.CASCADE)
    lig_chain = models.CharField(max_length=1)
    residues = models.JSONField(encoder=DjangoJSONEncoder)
    xtalform_site_id = models.TextField(
        null=False, help_text="xtalform site id from YAML"
    )
    # unused, remove when certain
    xtalform_site_num = models.TextField(
        null=True, help_text="alphabetic xtalform site id (enumerated on creation)"
    )

    objects = models.Manager()
    filter_manager = XtalformSiteDataManager()

    def __str__(self) -> str:
        return f"{self.xtalform_site_id}"

    def __repr__(self) -> str:
        return "<XtalformSite %r %r %r>" % (
            self.id,
            self.xtalform_site_id,
            self.xtalform,
        )


class CanonSiteConf(Versionable, models.Model):
    canon_site = models.ForeignKey(CanonSite, on_delete=models.CASCADE)
    # TODO: name not present in metadata atm
    name = models.TextField(null=True)
    ref_site_observation = models.OneToOneField(
        "SiteObservation", null=True, on_delete=models.CASCADE
    )
    residues = models.JSONField(encoder=DjangoJSONEncoder)

    objects = models.Manager()
    filter_manager = CanonSiteConfDataManager()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CanonSiteConf %r %r %r>" % (self.id, self.name, self.canon_site)


class Pose(models.Model):
    canon_site = models.ForeignKey(CanonSite, on_delete=models.CASCADE)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)
    main_site_observation = models.OneToOneField(
        "SiteObservation",
        null=True,
        on_delete=models.CASCADE,
        related_name="main_pose",
    )
    display_name = models.TextField(null=True)

    objects = models.Manager()
    filter_manager = PoseDataManager()

    def __str__(self) -> str:
        return f"{self.main_site_observation.code}"

    def __repr__(self) -> str:
        return "<Pose %r %r %r>" % (
            self.id,
            self.main_site_observation.code,
            self.display_name,
        )


class SiteObservation(Versionable, models.Model):
    SHORT_UUID_LENGTH: int = 4

    code = models.TextField(null=True)
    longcode = models.TextField(null=True)
    experiment = models.ForeignKey(Experiment, null=True, on_delete=models.CASCADE)
    cmpd = models.ForeignKey(Compound, null=True, on_delete=models.CASCADE)
    xtalform_site = models.ForeignKey(XtalformSite, null=True, on_delete=models.CASCADE)
    canon_site_conf = models.ForeignKey(
        CanonSiteConf, null=True, on_delete=models.CASCADE
    )
    pose = models.ForeignKey(
        Pose,
        on_delete=models.SET_NULL,
        null=True,
        related_name="site_observations",
    )
    bound_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    apo_solv_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    apo_desolv_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    apo_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    sigmaa_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    diff_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    event_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    artefacts_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    pdb_header_file = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    smiles = models.TextField()
    # rdkit representation of smiles field for structure-based
    # search. Internally rdkit mol type
    smiles_mol = models.TextField(editable=False, null=True)
    seq_id = models.IntegerField(null=True)
    chain_id = models.CharField(max_length=1, null=True)
    ligand_mol = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    ligand_smiles = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    ligand_pdb = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    ligand_sdf = models.FileField(
        upload_to="target_loader_data/", null=True, max_length=255
    )
    altloc = models.CharField(default='0', blank=True, max_length=1)

    # name is like 'v1a'
    # name generated for the virtual observation at upload
    virtual_name = models.TextField(null=True)

    # Set from the _Name property of the underlying Molecule
    virtual_molecule_name = models.TextField(null=True, blank=True)

    # A four character string of non-confusing uppercase letters and
    # digits for easy reference. This is combined with the Target to
    # form the ComputedMolecule's name",
    virtual_identifier = ShortUUIDField(
        length=SHORT_UUID_LENGTH,
        alphabet="ACDEFGHJKLMNPRSTUVWXYZ345679",
        null=True,
        blank=True,
    )

    # An optional url linking to the reference for this molecule
    virtual_ref_url = models.TextField(null=True, blank=True)

    # An optional rationale for this molecule
    virtual_rationale = models.TextField(null=True, blank=True)

    # Link to user-uploaded pdb file
    # NB! only uploaded pdb, not link to experiment.pdb_info, like before
    virtual_pdb_info = models.FileField(
        upload_to="computed_set_data/",
        null=True,
        max_length=255,
    )
    virtual_ligand_mol = models.FileField(
        upload_to="computed_set_data/", null=True, max_length=255
    )
    virtual_ref_observation = models.ForeignKey(
        'self', null=True, on_delete=models.CASCADE
    )

    computed_observations = models.ManyToManyField('self', blank=True)

    computed_inspirations = models.ManyToManyField(
        'self',
        through='ComputedInspiration',
        through_fields=("site_observation", "computed_inspiration"),
        symmetrical=False,
        related_name='inspired_observation',
    )

    objects = models.Manager()
    # causes problems with trigger func and don't really need it in
    # history anyway
    history = HistoricalRecords(excluded_fields=['smiles_mol'])
    filter_manager = SiteObservationDataManager()

    def __str__(self) -> str:
        return f"{self.code}"

    def __repr__(self) -> str:
        return "<SiteObservation %r %r %r %r>" % (
            self.id,
            self.code,
            self.experiment,
            self.cmpd,
        )

    def get_ligand_mol_file(self):
        contents = ''
        if self.ligand_mol:
            path = Path(settings.MEDIA_ROOT).joinpath(self.ligand_mol.name)
            with contextlib.suppress(TypeError, FileNotFoundError):
                with open(path, "r", encoding="utf-8") as f:
                    contents = f.read()

        return contents

    def get_filename(self):
        """Basename for this observation's uploaded pdb in downloads.

        Mirrors the former ComputedMolecule.get_filename: strip the
        auto-assigned suffix from virtual_pdb_info, e.g.
        ``computed_set_data/A0486a#<hash>.pdb_<hash>`` -> ``A0486a.pdb``.
        Returns None if there is no uploaded pdb.
        """
        if not self.virtual_pdb_info:
            return None
        fname = Path(self.virtual_pdb_info.name).name
        # With a referenced observation the name is already clean; without
        # one it still carries the auto-assigned '#<hash>' suffix to strip.
        if self.virtual_ref_observation:
            return fname
        if fname.find('#') > 0:
            return f"{fname.split('#')[0]}.pdb"
        return fname


class SiteObservationQualityStatus(models.Model):
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
    status = models.ForeignKey(QualityStatusType, on_delete=models.CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    timestamp = models.DateTimeField(default=timezone.now)
    auto_assigned = models.BooleanField(default=False)
    main_status = models.BooleanField(default=False)
    comment = models.TextField()

    objects = models.Manager()
    filter_manager = SiteObservationQualityStatusDataManager()

    def save(self, *args, **kwargs):
        try:
            with transaction.atomic():
                if self.main_status:
                    # lock rows for update to avoid race conditions
                    existing_main_statuses = (
                        SiteObservationQualityStatus.objects.select_for_update()
                        .filter(
                            site_observation=self.site_observation, main_status=True
                        )
                        .exclude(id=self.id)
                    )
                    existing_main_statuses.update(main_status=False)

                super().save(*args, **kwargs)
        except IntegrityError as e:
            # for some reason there's still a main status for this
            # observation. This is most probably a temporary glitch
            raise ValueError(
                "Another main_status already exists for site_observation "
                + f"{self.site_observation.id}"
            ) from e

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["site_observation"],
                condition=models.Q(main_status=True),
                name="unique_main_status_per_site_observation",
            )
        ]


class AtomCoordinates(models.Model):
    """Store ligand atom coordinates"""

    site_observation = models.ForeignKey(
        SiteObservation,
        on_delete=models.CASCADE,
        related_name="atom_coordinates",
    )
    atom_number = models.SmallIntegerField(null=False)
    coords = HalfVectorField(dimensions=3)

    class Meta:
        indexes = [
            HnswIndex(
                name='pgvector_coord_index',
                fields=['coords'],
                m=16,
                ef_construction=64,
                opclasses=['halfvec_l2_ops'],
            ),
        ]


class CompoundIdentifierType(models.Model):
    name = models.TextField(primary_key=True)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CompoundIdentifierType %r>" % (self.name)


class CompoundIdentifier(models.Model):
    type = models.ForeignKey(
        CompoundIdentifierType, to_field='name', on_delete=models.CASCADE
    )
    compound = models.ForeignKey(
        Compound,
        on_delete=models.CASCADE,
        related_name="all_identifiers",
    )
    url = models.URLField(null=True)
    name = models.TextField(null=False)

    objects = models.Manager()
    filter_manager = CompoundIdentifierDataManager()

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "type",
                    "compound",
                    "name",
                ],
                name="unique_compoundidentifier",
            ),
        ]

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CompoundIdentifier %r %r %r>" % (self.id, self.name, self.type)


class ActivityPoint(models.Model):
    source = models.CharField(max_length=50, null=True, db_index=True)
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    cmpd_id = models.ForeignKey(Compound, on_delete=models.CASCADE)
    activity = models.FloatField(db_index=True, help_text="Measured log(10) activity")
    units = models.CharField(max_length=50, help_text="Units (e.g. uM or whatever)")
    confidence = models.IntegerField(null=True, db_index=True)
    internal_id = models.CharField(
        max_length=150, null=True, help_text="The ID of the compound for internal use"
    )
    operator = models.CharField(
        max_length=5, default="NA", help_text="An operator, like > < or ="
    )

    def __str__(self) -> str:
        return f"{self.source}"

    def __repr__(self) -> str:
        return "<ActivityPoint %r %r %r %r %r %r>" % (
            self.id,
            self.source,
            self.target_id,
            self.activity,
            self.cmpd_id,
            self.units,
        )

    class Meta:
        unique_together = ("target_id", "activity", "cmpd_id", "units")


class ActionType(models.Model):
    id = models.AutoField(primary_key=True)
    description = models.CharField(max_length=200, default='')
    active = models.BooleanField(default=False)
    activation_date = models.DateTimeField(default=timezone.now)

    def __str__(self) -> str:
        return f"{self.description}"

    def __repr__(self) -> str:
        return "<ActionType %r %r %r>" % (self.id, self.description, self.active)

    class Meta:
        db_table = 'viewer_actiontype'


# Start of Session Project
class SessionProject(models.Model):
    title = models.CharField(max_length=200)
    init_date = models.DateTimeField(default=timezone.now)
    description = models.CharField(
        max_length=255,
        default='',
        help_text='A short user-defined description for the project',
    )
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    project = models.ForeignKey(
        Project,
        null=True,
        on_delete=models.CASCADE,
        help_text='Foreign Key link to the relevant project'
        ' (optional for legacy reasons)',
    )
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    tags = models.TextField(
        default='[]',
        help_text='A comma-separated list of user-defined tags'
        ' used for searching and tagging projects',
    )

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<SessionProject %r %r %r %r>" % (
            self.id,
            self.title,
            self.target,
            self.project,
        )

    class Meta:
        db_table = 'viewer_sessionproject'


class SessionActions(models.Model):
    """Django model for holding the user actions related to a particular session_project.

    actions is a JSON field containing types of actions related to the session_project.
    The list elements are (at the time of writing): -

    {
        "action_type" : "1",
        "action_datetime" : "2020-09-30T13:44:00.000Z",
        "object_type" : "",
        "object_name" : "",
        "show" : "true",
        "save" : "false",
    }
    """

    id = models.AutoField(primary_key=True)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    session_project = models.ForeignKey(SessionProject, on_delete=models.CASCADE)
    last_update_date = models.DateTimeField(default=timezone.now)
    actions = models.JSONField(encoder=DjangoJSONEncoder)

    objects = models.Manager()
    filter_manager = SessionActionsDataManager()

    def __str__(self) -> str:
        return f"{self.author}"

    def __repr__(self) -> str:
        return "<SessionActions %r %r %r>" % (
            self.id,
            self.author,
            self.session_project,
        )

    class Meta:
        db_table = 'viewer_sessionactions'


class Snapshot(models.Model):
    """Information describing a static snapshot of a fragalysis page.
    Multiple snapshots make up a project.
    """

    INIT = "INIT"
    AUTO = "AUTO"
    MANUAL = "MANUAL"
    SNAPSHOT_TYPE = (
        (INIT, "INIT"),  # Initial snapshot generated by system
        (AUTO, 'AUTO'),  # Automatic generated by system
        (MANUAL, 'MANUAL'),  # Manual generated by user action
    )
    id = models.AutoField(primary_key=True)
    type = models.CharField(choices=SNAPSHOT_TYPE, default=INIT, max_length=8)
    title = models.CharField(max_length=255)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    description = models.CharField(max_length=255, default='')
    created = models.DateTimeField(default=timezone.now)
    data = models.TextField(
        help_text='JSON data that is passed from the front-end'
        ' describing what to load into the react components'
        ' to reproduce the session state'
    )
    session_project = models.ForeignKey(
        SessionProject, null=True, on_delete=models.CASCADE
    )
    parent = models.ForeignKey(
        'self',
        models.DO_NOTHING,
        blank=True,
        null=True,
        related_name='children',
        help_text='Another Snapshot instance describing the current Snapshot parent',
    )
    additional_info = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        help_text='Optional JSON field containing name/value pairs for future use',
    )
    # NB! this field is accessed from a different serializer/endpoint
    state = models.JSONField(encoder=DjangoJSONEncoder, null=True)

    objects = models.Manager()
    filter_manager = SnapshotDataManager()

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<Snapshot %r %r %r %r>" % (self.id, self.title, self.type, self.author)

    class Meta:
        managed = True
        db_table = 'viewer_snapshot'


class SnapshotScreenshot(models.Model):
    snapshot = models.ForeignKey(Snapshot, null=False, on_delete=models.CASCADE)
    screenshot = models.TextField(null=True)
    screenshot_type = models.IntegerField(null=True)


class SnapshotActions(models.Model):
    """User actions leading to a particular snapshot or idea.

    'actions' is a JSON field containing types of actions made leading to the snapshot.
    The list elements are (at the time of writing): -

    {
        "action_type" : "1",
        "action_datetime" : "2020-09-30T13:44:00.000Z",
        "object_type" : "",
        "object_name" : "",
        "show" : "true",
        "save" : "false",
    }
    """

    id = models.AutoField(primary_key=True)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    session_project = models.ForeignKey(
        SessionProject, null=True, on_delete=models.CASCADE
    )
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE)
    last_update_date = models.DateTimeField(default=timezone.now)
    actions = models.JSONField(encoder=DjangoJSONEncoder)

    objects = models.Manager()
    filter_manager = SnapshotActionsDataManager()

    def __str__(self) -> str:
        return f"{self.author}"

    def __repr__(self) -> str:
        return "<SnapshotActions %r %r %r %r>" % (
            self.id,
            self.author,
            self.session_project,
            self.snapshot,
        )

    class Meta:
        db_table = 'viewer_snapshotactions'


class DesignSet(models.Model):
    """Information about a design sets - sets of 2D compounds that have been designed but do
    not yet have a 3D structure - unused
    """

    LIB = "library"
    FUP = "follow-up"
    USR = "user-submitted"
    ENM = "enumerated"
    SET_TYPE = (
        (LIB, "library"),  # library - e.g. DSiPoised
        (FUP, 'follow-up'),  # follow-up - e.g. purchased compounds
        (USR, 'user-submitted'),  # user submitted - can be submitted by anyone
        (ENM, 'enumerated'),  # enumerated - e.g. similarity search or something
    )
    compounds = models.ManyToManyField(
        Compound, help_text="The compounds that are in the design set"
    )
    set_name = models.CharField(max_length=50)
    set_type = models.CharField(max_length=100, choices=SET_TYPE, default=USR)
    set_description = models.TextField(max_length=1000, blank=True, null=True)

    def __str__(self) -> str:
        return f"{self.set_name}"

    def __repr__(self) -> str:
        return "<DesignSet %r %r %r>" % (self.id, self.set_name, self.set_type)


class ComputedSetSubmitter(models.Model):
    name = models.TextField()
    email = models.TextField()
    institution = models.TextField(
        help_text="The institution or organizational affiliation"
        " of the compound set submitter",
    )
    generation_date = models.DateField(null=True)
    method = models.TextField(
        help_text="A name for the method that was used" " to produce the uploaded data",
    )

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<ComputedSetSubmitter %r %r %r>" % (self.id, self.name, self.email)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "email",
                    "method",
                ],
                name="unique_computedsetsubmitter_email_method",
            ),
        ]


class CSetKeys(models.Model):
    """Used for authentication when uploading computed sets -
    each user is given an upload key associated with their email address in the form
    of a uuid. This is entered on the computed set upload page to allow a user upload.
    """

    user = models.CharField(max_length=50, default='User', editable=False)
    uuid = models.UUIDField(
        default=uuid.uuid4,
        editable=False,
        primary_key=True,
        help_text="Unique key for the user",
    )

    def __str__(self) -> str:
        return f"{self.user}"

    def __repr__(self) -> str:
        return "<CSetKeys %r %r>" % (self.uuid, self.user)


# computed sets = sets of poses calculated computationally
class ComputedSet(models.Model):
    """Computed sets - sets of 3D poses of molecules calculated computationally
    and uploaded by a user
    """

    # if null, then resultupload
    name = models.TextField(null=True)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE)
    submitted_sdf = models.FileField(
        upload_to='computed_set_data/',
        help_text="The original SDF containing the ComputedSet",
    )
    written_sdf_filename = models.TextField(
        null=True,
        help_text="The written ComputedSet filename",
    )
    spec_version = models.FloatField(
        null=True, help_text="The version of the SDF file format specification"
    )
    method_url = models.TextField(
        null=True,
        help_text="A url linking to a write-up of the methodology used to create the"
        " computed set",
    )
    submitter = models.ForeignKey(
        ComputedSetSubmitter, null=True, on_delete=models.CASCADE
    )
    method = models.TextField(
        null=True,
        blank=True,
        help_text="The name of the algorithmic method used to generate the compounds (e.g. Fragmenstein)",
    )
    upload_date = models.DateField(
        null=True,
        blank=True,
        help_text="The date the set was uploaded",
    )
    md_ordinal = models.SmallIntegerField(
        null=True,
        blank=True,
        help_text="The ordinal distinguishing between uploads using the same method and date",
    )
    # RU: uploaded_by. this is confusing, submitter sounds like it
    # could be better, but it has it's own model
    owner_user = models.ForeignKey(
        User, on_delete=models.CASCADE, default=settings.ANONYMOUS_USER
    )
    upload_datetime = models.DateTimeField(
        null=True,
        blank=True,
        default=timezone.now,
    )
    site_observations = models.ManyToManyField(
        SiteObservation,
        through="ComputedSetSiteObservation",
        through_fields=("computed_set", "site_observation"),
        related_name="computed_set",
    )

    objects = models.Manager()
    filter_manager = ComputedSetDataManager()
    history = HistoricalRecords()

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "name",
                    "target",
                ],
                name="unique_computedsetname_target",
            ),
        ]

    def __str__(self) -> str:
        target_title: str = self.target.title if self.target else "None"
        return f"{self.name} {target_title}"

    def __repr__(self) -> str:
        return "<ComputedSet %r %r %r>" % (self.id, self.name, self.target)


class ComputedInspiration(models.Model):
    site_observation = models.ForeignKey(
        SiteObservation,
        on_delete=models.CASCADE,
        related_name='+',
    )
    computed_inspiration = models.ForeignKey(
        SiteObservation,
        on_delete=models.CASCADE,
        related_name='+',
    )
    computed_set = models.ForeignKey(
        ComputedSet,
        on_delete=models.CASCADE,
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "site_observation",
                    "computed_inspiration",
                    "computed_set",
                ],
                name="unique_inspirations_computed_set",
            ),
        ]


class ComputedSetSiteObservation(models.Model):
    computed_set = models.ForeignKey(ComputedSet, null=False, on_delete=models.CASCADE)
    site_observation = models.ForeignKey(
        SiteObservation, null=False, on_delete=models.CASCADE
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "computed_set",
                    "site_observation",
                ],
                name="unique_computedsetsiteobservation",
            ),
        ]


class SiteObservationComputedSiteObservation(models.Model):
    """Store alignment matches between experimental and computed observations.

    On upload, look for uploaded ComputedSets and observations
    (formerly ComputedMolecules, RHS compounds), compare the
    alignments and store the matches found along with the RMSD.

    """

    site_observation = models.ForeignKey(
        SiteObservation,
        null=False,
        on_delete=models.CASCADE,
        related_name="lhs_site_observations",
    )
    computed_site_observation = models.ForeignKey(
        SiteObservation,
        null=False,
        on_delete=models.CASCADE,
        related_name="rhs_site_observations",
    )
    rmsd = models.FloatField(null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "site_observation",
                    "computed_site_observation",
                ],
                name="unique_siteobservation_computedsiteobservation",
            ),
        ]


class File(models.Model):
    file = models.FileField(blank=False)

    def __str__(self):
        return self.file.name

    def __repr__(self) -> str:
        return "<File %r %r>" % (self.id, self.file.name)


class DiscourseCategory(models.Model):
    """Discourse Subcategory references for Fragalysis - initially Targets"""

    category_name = models.CharField(
        max_length=200,
        unique=True,
        help_text="The name of the (sub)category within Discourse. It must be unique",
    )
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    discourse_category_id = models.IntegerField(
        help_text="The Discourse categoryID." " Returned when the category was created"
    )

    def __str__(self):
        return self.author.username

    def __repr__(self) -> str:
        return "<DiscourseCategory %r %r %r>" % (
            self.id,
            self.author,
            self.category_name,
        )

    class Meta:
        db_table = 'viewer_discoursecategory'


class DiscourseTopic(models.Model):
    """Discourse Topic references for Fragalysis - initially Targets"""

    topic_title = models.CharField(
        max_length=200,
        unique=True,
        validators=[
            MinLengthValidator(
                15, 'Discourse Topic Title must be longer than 15 characters'
            )
        ],
        help_text="The title of the (sub)category within Discourse."
        " It must be unique within Discourse",
    )
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    discourse_topic_id = models.IntegerField()

    def __str__(self):
        author: str = self.author.username if self.author else "None"
        return f"{author} '{self.topic_title}'"

    def __repr__(self) -> str:
        return "<DiscourseTopic %r %r %r>" % (self.id, self.author, self.topic_title)

    class Meta:
        db_table = 'viewer_discoursetopic'


class DownloadLinks(models.Model):
    """Searches made with the download_structures api."""

    # Stores the basename of the download file (e.g. "TARGET.zip"). The
    # absolute path is reconstructed via get_file_url() using MEDIA_ROOT,
    # the "downloads" subdir and task_id. Not unique because two records
    # for different tasks can produce the same filename. Indexed via
    # Meta.indexes below (looked up by basename in download dedup).
    file_url = models.TextField(null=True)
    task_id = models.TextField(
        null=True,
        help_text="The task ID assigned to this download (if a Task is launched)",
    )
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    target = models.ForeignKey(
        Target, null=True, on_delete=models.CASCADE, db_index=True
    )
    # list of sorted observation shortcodes. Changed in 1982, used to
    # be JSONfield but this didn't work with .get(). The same fate may
    # wait for the other json fields, they're only needed for
    # comparison and the content's isn't really
    # Update: changed again in 2142 to contain ids instead of names
    proteins = models.TextField(null=True)
    protein_params = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        help_text="Contains a sorted list of parameters used to create the zip file",
    )
    other_params = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        help_text="Contains a sorted list of parameters used to create the zip file",
    )
    static_link = models.BooleanField(
        default=False, help_text="This preserves the proteins from the previous search"
    )
    # TODO - zip_contents is no longer Used (A.Christie 2024-01-19)
    zip_contents = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        help_text="For static files, this field contains the contents of the zip so that"
        " it can be reconstructed with the same file-links that it had previously."
        " For dynamic files, the zip is reconstructed from the search",
    )
    create_date = models.DateTimeField()
    expired_date = models.DateTimeField(
        null=True,
        help_text="Set when the download has passed its keep until date."
        " Users should not use records when this is set.",
    )
    keep_zip_until = models.DateTimeField(
        db_index=True,
        null=True,
        help_text="The datetime when the tag was created"
        " plus the retention time"
        " (1 hour at the time of writing)",
    )
    deleted = models.BooleanField(
        null=True,
        help_text="Set when the download file has been removed from the filesystem.",
    )
    expiry_reason = models.TextField(
        null=True,
        help_text="Why the download was expired (e.g. the task was lost)."
        " Shown to users querying a failed download.",
    )
    # TODO - zip_file is no longer Used (A.Christie 2024-01-19)
    zip_file = models.BooleanField(default=False)
    original_search = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    request_ip = models.TextField(null=True)
    request_location = models.TextField(null=True)

    def __str__(self):
        return str(self.file_url)

    def __repr__(self) -> str:
        return "<DownloadLinks %r %r %r>" % (
            self.file_url,
            self.user,
            self.target,
        )

    def get_file_url(self):
        """Reconstruct the absolute filesystem path of the download file.

        Returns None when either file_url (the basename) or task_id (the
        per-download directory name) is missing — i.e. the download has
        not been built or has been hard-erased.
        """
        if not self.file_url or not self.task_id:
            return None
        return os.path.join(
            settings.MEDIA_ROOT, 'downloads', self.task_id, self.file_url
        )

    class Meta:
        db_table = 'viewer_downloadlinks'
        indexes = [
            models.Index(
                fields=['file_url'],
                name='downloadlinks_file_url_idx',
            ),
        ]


class TagCategory(models.Model):
    category = models.CharField(
        max_length=50, unique=True, help_text="The name of the tag category"
    )
    colour = models.CharField(
        max_length=20, null=True, help_text="Expected to be an RGB string"
    )
    description = models.CharField(max_length=200, null=True)

    def __str__(self):
        return str(self.category)

    def __repr__(self) -> str:
        return "<TagCategory %r %r>" % (self.id, self.category)

    class Meta:
        db_table = 'viewer_tagcategory'


class Tag(models.Model):
    tag = models.TextField(help_text="The (unique) name of the tag")
    short_tag = models.TextField(
        null=True,
        help_text="The generated shorter version of tag (without target name)",
    )
    tag_prefix = models.TextField(
        null=True, help_text="Tag prefix for auto-generated tags"
    )
    upload_name = models.TextField(null=True)
    category = models.ForeignKey(TagCategory, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    create_date = models.DateTimeField(default=timezone.now)
    colour = models.CharField(
        max_length=20, null=True, help_text="Expected to be an RGB string"
    )
    discourse_url = models.TextField(max_length=1000, null=True)
    help_text = models.TextField(null=True)
    meta_category = models.TextField(null=True)
    additional_info = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        help_text="Optional JSON field containing name/value pairs for future use",
    )
    hidden = models.BooleanField(default=False)

    def __str__(self) -> str:
        return f"{self.tag}"

    def __repr__(self) -> str:
        return "<Tag %r %r %r %r %r %r>" % (
            self.id,
            self.tag,
            self.upload_name,
            self.category,
            self.target,
            self.user,
        )

    class Meta:
        abstract = True
        unique_together = (
            'upload_name',
            'target',
        )


class SiteObservationTag(Tag):
    site_observations = models.ManyToManyField(
        SiteObservation,
        through="SiteObvsSiteObservationTag",
        through_fields=("site_obvs_tag", "site_observation"),
    )
    mol_group = models.ForeignKey(
        "scoring.SiteObservationGroup", null=True, blank=True, on_delete=models.SET_NULL
    )

    def __str__(self) -> str:
        return f"{self.id}"

    def __repr__(self) -> str:
        return "<SiteObservationTag %r %r>" % (self.id, self.site_observations)


class SiteObvsSiteObservationTag(models.Model):
    site_obvs_tag = models.ForeignKey(
        SiteObservationTag, null=False, on_delete=models.CASCADE
    )
    site_observation = models.ForeignKey(
        SiteObservation, null=False, on_delete=models.CASCADE
    )

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "site_observation",
                    "site_obvs_tag",
                ],
                name="unique_siteobservationtagcontents",
            ),
        ]


class SessionProjectTag(Tag):
    """Data for SessionProjectTag(s) inherited from Tag."""

    session_projects = models.ManyToManyField(SessionProject)

    def __str__(self) -> str:
        return f"{self.id}"

    def __repr__(self) -> str:
        return "<SessionProjectTag %r>" % self.id


class JobFileTransfer(models.Model):
    PENDING = "PENDING"
    STARTED = "STARTED"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"
    RETRY = "RETRY"
    REVOKED = "REVOKED"
    STATUS = (
        (PENDING, 'PENDING'),  # Initial state when queued
        (STARTED, 'STARTED'),  # File transfer started
        (SUCCESS, 'SUCCESS'),  # File transfer finished successfully
        (FAILURE, 'FAILURE'),  # File transfer failed
        (RETRY, 'RETRY'),
        (REVOKED, 'REVOKED'),
    )
    id = models.AutoField(primary_key=True)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE)
    target = models.ForeignKey(
        Target, null=True, on_delete=models.CASCADE, db_index=True
    )
    squonk_project = models.CharField(max_length=200, null=True)
    sub_path = ShortUUIDField(
        length=4, alphabet="abcdefghijklmnopqrstuvwxyz", null=True
    )
    # A list of files (paths and files relative to MEDIA_ROOT)...
    proteins = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    # A list of files (paths and files relative to MEDIA_ROOT)...
    compounds = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    transfer_task_id = models.CharField(null=True, max_length=50)
    transfer_status = models.CharField(choices=STATUS, default=PENDING, max_length=7)
    transfer_progress = models.DecimalField(
        null=True,
        max_digits=5,
        decimal_places=2,
        help_text="Intended to be used as an indication of progress (0 to 100%)",
    )
    transfer_datetime = models.DateTimeField(
        null=True, help_text="The datetime the transfer was completed"
    )

    def __str__(self) -> str:
        return f"{self.user}"

    def __repr__(self) -> str:
        return "<JobFileTransfer %r %r %r %r %r>" % (
            self.id,
            self.user,
            self.snapshot,
            self.target,
            self.squonk_project,
        )

    class Meta:
        db_table = 'viewer_jobfiletransfer'


class JobRequest(models.Model):
    PENDING = "PENDING"
    STARTED = "STARTED"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"
    RETRY = "RETRY"
    REVOKED = "REVOKED"
    SQUONK_STATUS = (
        (PENDING, 'PENDING'),  # Initial state when job queued
        (STARTED, 'STARTED'),  # Job started in squonk (updated by squonk)
        (
            SUCCESS,
            'SUCCESS',
        ),  # Job completed successfully in squonk (updated by squonk)
        (FAILURE, 'FAILURE'),  # Job failed in squonk (updated by squonk)
        (RETRY, 'RETRY'),  # Job status in squonk (updated by squonk)
        (REVOKED, 'REVOKED'),  # Job status in squonk (updated by squonk)
    )
    UPLOAD_STATUS = (
        (PENDING, 'PENDING'),  # Initial state when upload queued
        (STARTED, 'STARTED'),  # Upload job started
        (SUCCESS, 'SUCCESS'),  # Upload job successful
        (FAILURE, 'FAILURE'),  # Upload job failed
        (RETRY, 'RETRY'),
        (REVOKED, 'REVOKED'),
    )
    id = models.AutoField(primary_key=True)
    squonk_job_name = models.CharField(max_length=200, null=True)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    snapshot = models.ForeignKey(
        Snapshot,
        on_delete=models.CASCADE,
        help_text="The snapshot the file transfer is part of",
    )
    target = models.ForeignKey(
        Target, null=True, on_delete=models.CASCADE, db_index=True
    )
    project = models.ForeignKey(Project, null=True, on_delete=models.CASCADE)
    squonk_project = models.CharField(
        max_length=200,
        null=True,
        help_text="The name of a project that has been created in Squonk"
        " that the files will be transferred to",
    )
    squonk_job_spec = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    job_start_datetime = models.DateTimeField(null=True)
    job_finish_datetime = models.DateTimeField(
        null=True,
        help_text="The datetime when the Squonk Job has finished,"
        " populated by information in the Squonk callback."
        " If this is not set you can assume the JOb is still"
        " running. When it is set the job_status filed will be"
        " updated (to SUCCESS or FAILURE). If automatic upload"
        " follows an upload_task_id wil be set and you can"
        " monitor upload_status for a status of the upload",
    )
    job_status = models.CharField(
        choices=SQUONK_STATUS,
        default=PENDING,
        max_length=7,
        help_text="The status of the Squonk job, e.g. 'PENDING'"
        " Will be modified by Squonk through the callback URL",
    )
    job_status_datetime = models.DateTimeField(null=True)
    squonk_job_info = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        help_text="squonk_job_info is a copy of the response from DmApi.start_job_instance()."
        " It's an instance of a DmApiRv object (a namedtuple)"
        " that contains a 'success' (boolean) and 'msg' (the DmApi response's resp.json())."
        " For us this will contain a 'task_id', 'instance_id' and 'callback_token'."
        " The content will be a list with index '0' that's the value of the DmApiRv"
        " 'success' variable and, at index '1', the original response message json()."
        " The Job callback token will be squonk_job_info[1]['callback_token']",
    )
    squonk_url_ext = models.CharField(
        max_length=200,
        null=True,
        help_text="a Squonk UI URL to obtain information about the running instance."
        " It's essentially the Squonk URL with the instance ID appended.",
    )
    code = models.UUIDField(
        default=uuid.uuid4,
        editable=False,
        unique=True,
        help_text="A UUID generated by Fragalysis and passed to Squonk as part of a callback URL."
        " This value is used to uniquely identify the HJob in Squonk and is passed back"
        " by squonk to provide context in calls to the JobCallBackView",
    )
    upload_task_id = models.CharField(
        null=True,
        max_length=50,
        help_text="Celery task ID for results upload task."
        " Set when the Job completes and an automated upload follows",
    )
    upload_status = models.CharField(
        choices=UPLOAD_STATUS,
        default=PENDING,
        max_length=7,
        null=True,
        help_text="Status of upload task",
    )
    computed_set = models.ForeignKey(ComputedSet, on_delete=models.CASCADE, null=True)

    def __str__(self) -> str:
        return f"{self.user}"

    def __repr__(self) -> str:
        return "<JobRequest %r %r %r %r %r %r>" % (
            self.id,
            self.user,
            self.squonk_job_name,
            self.snapshot,
            self.target,
            self.squonk_project,
        )

    class Meta:
        db_table = 'viewer_jobrequest'


class JobOverride(models.Model):
    override = models.JSONField(encoder=DjangoJSONEncoder)
    author = models.ForeignKey(
        User,
        null=True,
        on_delete=models.SET_NULL,
        help_text="The user that uploaded the override",
    )

    def __str__(self) -> str:
        return f"{self.author}"

    def __repr__(self) -> str:
        return "<JobOverride %r %r>" % (self.id, self.author)

    class Meta:
        db_table = 'viewer_joboverride'


class Squonk2Org(models.Model):
    """Squonk2 Organisations (UUIDs) and the Account Servers
    they belong to. Managed by the Squonk2Agent class and only one entry expected.
    """

    uuid = models.TextField(
        max_length=40,
        help_text="A Squonk2 Account Server (AS) Organisation UUID."
        " A fixed length string consisting of 'org-' followed by a uuid4 value,"
        " e.g. 'org-54260047-183b-42e8-9658-385a1e1bd236'",
    )
    name = models.TextField(
        max_length=80,
        help_text="The name of the Squonk2 Organisation UUID (obtained form the AS)",
    )
    as_url = models.URLField()
    as_version = models.TextField()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Squonk2Org %r %r %r>" % (self.id, self.name, self.uuid)


class Squonk2Unit(models.Model):
    """Squonk2 Unit (UUIDs). Managed by the Squonk2Agent class."""

    uuid = models.TextField(
        max_length=41,
        help_text="A Squonk2 Account Server (AS) Unit UUID."
        " A fixed length string consisting of 'unit-' followed by a uuid4 value,"
        " e.g. 'unit-54260047-183b-42e8-9658-385a1e1bd236'",
    )
    name = models.TextField(
        help_text="The name used to create the Squonk2 Unit UUID"
        " This is not limited by the actual name length imposed by the DM"
    )
    organisation = models.ForeignKey(Squonk2Org, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Squonk2Unit %r %r %r>" % (self.id, self.name, self.uuid)


class Squonk2Project(models.Model):
    """Squonk2 Project and Product (UUIDs). Managed by the Squonk2Agent class."""

    uuid = models.TextField(max_length=44)
    name = models.TextField(
        help_text="The name of the Squonk2 Project UUID (obtained form the AS)"
    )
    product_uuid = models.TextField(
        max_length=44,
        help_text="A Squonk2 Account Server (AS) Product UUID."
        " A fixed length string consisting of 'product-' followed by a uuid4 value,"
        " e.g. 'product-54260047-183b-42e8-9658-385a1e1bd236'",
    )
    unit = models.ForeignKey(Squonk2Unit, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Squonk2Project %r %r %r %r %r>" % (
            self.id,
            self.name,
            self.uuid,
            self.product_uuid,
            self.unit,
        )


class ResultValueDataType(models.Model):
    data_type = models.TextField(primary_key=True)
    db_type = models.TextField()


class ResultValueModifier(models.Model):
    modifier = models.TextField(primary_key=True)
    django_operator = models.TextField()


class ResultUpload(models.Model):
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE)
    upload_file = models.FileField(upload_to='assay_data/')
    upload_date = models.DateTimeField(
        null=True,
        blank=True,
        default=timezone.now,
    )
    uploaded_by = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        default=settings.ANONYMOUS_USER,
    )

    objects = models.Manager()
    filter_manager = ResultUploadDataManager()


class ResultProperty(models.Model):
    """Assay data property name"""

    result_property = models.TextField()
    unit = models.TextField(null=True)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE)
    visible = models.BooleanField(default=True, null=False)
    order = models.PositiveSmallIntegerField(null=False, default=0, blank=True)
    data_type = models.ForeignKey(ResultValueDataType, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.result_property}"

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "result_property",
                    "unit",
                    "target",
                ],
                name="unique_result_property_target_unit",
            ),
        ]


class Result(models.Model):
    """Model to store assay results data"""

    raw_value = models.TextField(null=True)
    float_value = models.FloatField(null=True)
    int_value = models.IntegerField(null=True)
    link_value = models.TextField(null=True)
    numeric_modifier = models.ForeignKey(
        ResultValueModifier,
        on_delete=models.CASCADE,
        null=True,
    )
    text_value = models.TextField(null=True)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE, null=True)
    site_observation = models.ForeignKey(
        SiteObservation,
        on_delete=models.CASCADE,
        null=True,
    )
    experiment = models.ForeignKey(Experiment, null=True, on_delete=models.CASCADE)
    result_upload = models.ForeignKey(
        ResultUpload,
        on_delete=models.CASCADE,
        null=True,
    )
    # replacing result_upload with computed_set
    computed_set = models.ForeignKey(
        ComputedSet,
        on_delete=models.CASCADE,
        null=True,
    )
    parsing_error = models.BooleanField(default=False, null=False)
    result_property = models.ForeignKey(ResultProperty, on_delete=models.CASCADE)

    objects = models.Manager()
    filter_manager = AssayResultDataManager()

    def __str__(self) -> str:
        return f"{self.id}: {self.raw_value} {self.result_property.data_type}"


class PlotDataIdentifierType(models.Model):
    identifier = models.TextField(primary_key=True)


class PlotData(models.Model):
    """Store uploaded plotly plot data"""

    author = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        default=settings.ANONYMOUS_USER,
    )
    title = models.TextField()
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    identifier = models.ForeignKey(PlotDataIdentifierType, on_delete=models.CASCADE)
    upload_time = models.DateTimeField(
        blank=True,
        default=timezone.now,
    )
    plotly_data = models.JSONField(
        encoder=DjangoJSONEncoder,
        null=True,
        blank=True,
    )
    notebook_path = models.TextField(null=True)
    squonk_project_id = models.TextField(null=True)
