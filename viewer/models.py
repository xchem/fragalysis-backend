import uuid
from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone
from django.core.serializers.json import DjangoJSONEncoder
from django.core.validators import MinLengthValidator
from django.conf import settings

from django.contrib.postgres.fields import ArrayField

from simple_history.models import HistoricalRecords

from viewer.target_set_config import get_mol_choices, get_prot_choices

class Project(models.Model):
    title = models.CharField(max_length=200, unique=True)
    init_date = models.DateTimeField(auto_now_add=True)
    user_id = models.ManyToManyField(User)
    open_to_public = models.BooleanField(default=False)

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<Project %r %r %r>" % (self.id, self.title, self.open_to_public)


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
        (REVOKED, 'REVOKED')
    )

    title = models.CharField(unique=True, max_length=200, help_text="A title, i.e. Mpro")
    init_date = models.DateTimeField(auto_now_add=True)
    project_id = models.ManyToManyField(Project)
    uniprot_id = models.CharField(max_length=100, null=True,
                                  help_text="The uniprot ID id for the target. A unique key")
    metadata = models.FileField(upload_to="metadata/", null=True, max_length=255,
                                help_text='Optional file upload defining metadata about the target.'
                                          ' Can be used to add custom site labels')
    zip_archive = models.FileField(upload_to="archive/", null=True, max_length=255,
                                   help_text='Link to zip file created from targets uploaded with the loader')
    default_squonk_project = models.CharField(max_length=200, null=True)
    upload_task_id = models.CharField(null=True, max_length=50,
                                      help_text='The Task ID of upload Celery task)')
    upload_status = models.CharField(choices=STATUS, null=True, max_length=7,
                                     help_text='Identifies the status of the upload.'
                                               ' Will only be updated at the end of the process')
    upload_progress = models.DecimalField(null=True, max_digits=5, decimal_places=2,
                                          help_text='Intended to be used as an indication of progress (0 to 100%)')
    upload_datetime = models.DateTimeField(null=True,
                                           help_text='The datetime the upload was completed')

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<Target %r %r %r>" % (self.id, self.title, self.project_id)


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
    commit_datetime = models.DateTimeField(help_text="The UTC datetime the upload was committed")
    committer = models.ForeignKey(User, on_delete=models.CASCADE,
                                  help_text="The user committing the original file."
                                            " This user may not be the author of the file")
    complete_datetime = models.DateTimeField(null=True,
                                             help_text="The UTC datetime the upload finished."
                                                       " It can be considered a success"
                                                       " if the status is LOADED")
    task_id = models.CharField(null=True, max_length=50,
                               help_text="Celery task ID responsible for the upload")
    status = models.CharField(choices=STATUS, default=LOADING, max_length=7)
    message = models.TextField(null=True, blank=True,
                               help_text="Any message associated with the upload."
                                         " Typically set when status is FAILURE")

    def __str__(self) -> str:
        return f"{self.project}"

    def __repr__(self) -> str:
        return "<ExperimentUpload %r %r %r>" % (self.id, self.project, self.target)


class Experiment(models.Model):
    experiment_upload = models.ForeignKey(ExperimentUpload, on_delete=models.CASCADE)
    code = models.TextField(null=True)
    status = models.IntegerField(null=True)
    version = models.IntegerField(null=True)
    pdb_info = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    mtz_info = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    cif_info = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    event_map_info = ArrayField(models.FileField(), null=True)
    type = models.PositiveSmallIntegerField(null=True)
    pdb_sha256 = models.TextField(null=True)
    compounds = models.ManyToManyField("Compound")

    def __str__(self) -> str:
        return f"{self.code}"

    def __repr__(self) -> str:
        return "<Experiment %r %r %r>" % (self.id, self.code, self.experiment_upload)


# TODO: delete
class Protein(models.Model):
    """Information about a protein. A protein is a protein structure which has a unique set of
    3D coordinates, rather than a target, which is a set of protein objects of the same protein.
    A Molecule object is also linked to a protein, so that a complete structure is
    comprised of the molecule and protein in separate parts in fragalysis.

    protein type is from a pre-defined list and determined by file extension on upload
    (defined in loader.config.get_prot_choices): -

    prot_choices = (
        (APO, "Apo", "_apo.pdb", "APO"),
        (STRIPPED, "Stripped", "_no_buffer_altlocs.pdb", "STRIPPED"),
        (TLEAPED, "Tleaped", "_tleap.pdb", "TLEAP"),
        (CHUNKED, "Chunked", "_chunk.pdb", "CHUNK"),
        (BOUND, "Bound", '_bound.pdb', "BOUND")
    )
    """
    code = models.CharField(max_length=50, db_index=True, help_text='Code for this protein (e.g. NUDT5A-x0001_1A)')
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    experiment = models.ForeignKey(Experiment, null=True, on_delete=models.CASCADE)
    apo_holo = models.BooleanField(null=True, help_text='0 for apo (ligand removed), 1 for holo (ligand in tact)')
    # Set the groups types
    prot_choices, default_prot = get_prot_choices()
    prot_type = models.CharField(
        choices=prot_choices, default=default_prot, max_length=2
    )
    pdb_info = models.FileField(upload_to="pdbs/", null=True, max_length=255,
                                help_text='File link to apo pdb structure - pdb file with ligand removed')
    bound_info = models.FileField(upload_to="bound/", null=True, max_length=255,
                                  help_text='File link to uploaded bound pdb file (optional)')
    cif_info = models.FileField(upload_to="cifs/", null=True, max_length=255,
                                help_text='File link to uploaded cif file (optional)')
    mtz_info = models.FileField(upload_to="mtzs/", null=True, max_length=255,
                                help_text='File link to uploaded mtz file (optional)')
    map_info = models.FileField(upload_to="maps/", null=True, max_length=255,
                                help_text='File link to uploaded map file (optional)')
    sigmaa_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    diff_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    event_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    trans_matrix_info = models.FileField(upload_to="trans/", null=True, max_length=255,
                                         help_text='File link to uploaded transformation matrix file (optional)')
    pdb_header_info = models.FileField(upload_to="pdbs/", null=True, max_length=255,
                                       help_text='File link to uploaded _header.pdb file (optional)')
    apo_desolve_info = models.FileField(upload_to="pdbs/", null=True, max_length=255,
                                        help_text='File link to uploaded _apo-desolv.pdb file (optional)')
    aligned = models.BooleanField(null=True)
    aligned_to = models.ForeignKey("self", null=True, on_delete=models.CASCADE,
                                   help_text='Foreign key to another instance of a protein'
                                             ' to which this protein is aligned (optional)')
    has_eds = models.BooleanField(null=True)

    class Meta:
        unique_together = ("code", "target_id", "prot_type")

    def __str__(self):
        return self.code

    def __repr__(self):
        return "<Protein %s %s %s %s>" % (self.id, self.code, self.target_id, self.prot_type)


class Compound(models.Model):
    """Information about a compound, which is a unique 2D molecule
    """
    inchi = models.TextField(unique=False, db_index=True)
    smiles = models.CharField(max_length=255, db_index=True)
    current_identifier = models.CharField(max_length=255, db_index=True, blank=True, null=True,
                                          help_text='The identifier for this compound that is used in Fragalysis to'
                                                    ' represent its 3D molecule (optional)')
    all_identifiers = models.TextField(blank=True, null=True,
                                       help_text='A comma separated list of all identifiers that have been used in'
                                                 ' the past to represent this 2D compound')
    project_id = models.ManyToManyField(Project)
    mol_log_p = models.FloatField(help_text='Computed LogP value')
    mol_wt = models.FloatField(help_text='Computed molecular weight (Da)')
    tpsa = models.FloatField(help_text='Computed Topological Polar Surface Area')
    heavy_atom_count = models.IntegerField(help_text='Computed heavy (non-hydrogen) atom count')
    heavy_atom_mol_wt = models.FloatField(help_text='Computed molecular weight of all heavy (non-hydrogen) atoms')
    nhoh_count = models.IntegerField(help_text='Computed number of hydroxylamine groups')
    no_count = models.IntegerField(help_text='Computed number of nitroso groups')
    num_h_acceptors = models.IntegerField(help_text='Computed number of hydrogen-bond acceptor groups')
    num_h_donors = models.IntegerField(help_text='Computed number of hydrogen-bond donor groups')
    num_het_atoms = models.IntegerField(help_text='Computed number of heterogeneous atoms')
    num_rot_bonds = models.IntegerField(help_text='Computed number of rotatable bonds')
    num_val_electrons = models.IntegerField(help_text='Computed number of valence electrons')
    ring_count = models.IntegerField(help_text='Computed number of rings in the molecule')
    inspirations = models.ManyToManyField("Molecule", blank=True,
                                          help_text='Foreign key link to any number of 3D Molecules that inspired'
                                                    ' the design of this compound')
    description = models.TextField(blank=True, null=True)
    comments = models.TextField(blank=True, null=True)

    def __str__(self) -> str:
        return f"{self.smiles}"

    def __repr__(self) -> str:
        return "<Compound %r %r %r>" % (self.id, self.smiles, self.inchi)


# TODO: delete
class Molecule(models.Model):
    """Information about a Molecule. A molecule is linked to a compound, and represents the
    unique 3D coordinates for that molecule. Note a compound can be linked to multiple molecules,
    each with a different 3D structure, but the same 2D structure

    mol_type is from a pre-defined list and determined by file extension on upload
    (defined in loader.config.get_mol_choices):

    mol_choices = (
        (PROASIS, "Proasis molecule", ".mol", "MOL"),
        (SDF, "Sdf molecule", ".sdf", "SDF"),
        (HYDROGEN, "Hydrogens added ", "_h.mol", "H_MOL"),
        (
            HYDROGEN_AM1_CHARGES,
            "Mol2 format with Hydrogens and AM1 BCC",
            ".mol2",
            "MOL2",
        ),
        )
    )
    """
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    lig_id = models.CharField(max_length=5, null=True,
                              help_text="Ligand ID from the pdb file it originated from (e.g. LIG 1 A)")
    chain_id = models.CharField(max_length=1, null=True,
                                help_text="Chain ID from the pdb file it originated from (e.g. A)")
    mol_choices, default_mol = get_mol_choices()
    mol_type = models.CharField(choices=mol_choices, default=default_mol, max_length=2)
    sdf_info = models.TextField(null=True, help_text="The 3D coordinates for the molecule in MDL (mol file) format."
                                                     " Taken directly from the uploaded file")
    rscc = models.FloatField(null=True, help_text="The RSCC score of the molecule 3D coordinates vs. PANDDA event map")
    occupancy = models.FloatField(null=True, help_text="The occupancy (electron density) of the molecule")
    x_com = models.FloatField(null=True, help_text="x-coordinate for centre of mass")
    y_com = models.FloatField(null=True, help_text="y-coordinate for centre of mass")
    z_com = models.FloatField(null=True, help_text="z-coordinate for centre of mass")
    rmsd = models.FloatField(null=True)
    prot_id = models.ForeignKey(Protein, on_delete=models.CASCADE,
                                help_text="Associated protein (apo) that this ligand was pulled from")
    cmpd_id = models.ForeignKey(Compound, on_delete=models.CASCADE,
                                help_text="Associated 2D compound")
    sdf_file = models.FileField(upload_to="sdfs/", null=True, max_length=255,
                                help_text="File link to uploaded sdf file")
    # Tracks the changes made to an instance of this model over time
    history = HistoricalRecords()

    def __str__(self) -> str:
        return f"{self.smiles}"

    def __repr__(self) -> str:
        return "<Molecule %r %r %r %r %r>" % (self.id, self.smiles, self.mol_type, self.prot_id, self.cmpd_id)

    class Meta:
        unique_together = ("prot_id", "cmpd_id", "mol_type")



class QuatAssembly(models.Model):
    chains = models.TextField()
    name = models.TextField()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<QuatAssembly %r %r %r>" % (self.id, self.name, self.chains)


class Xtalform(models.Model):
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)
    # TODO: add fk to ref experiment.
    # DONE? this is what experiment is, do I need to change it?
    name = models.TextField(null=True)
    quat_assembly = models.ForeignKey(QuatAssembly, on_delete=models.CASCADE, null=True)
    space_group = models.TextField(null=True)
    unit_cell_info = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    xtalform_id = models.IntegerField(null=True, help_text="xtalform id from YAML")

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Xtalform %r %r %r>" % (self.id, self.name, self.experiment)


class CompoundIdentifierType(models.Model):
    NAME_LENGTH = 20
    name = models.CharField(max_length=NAME_LENGTH)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CompoundIdentifierType %r %r>" % (self.id, self.name)


class CompoundIdentifier(models.Model):
    NAME_LENGTH = 40
    URL_LENGTH = 200
    type = models.ForeignKey(CompoundIdentifierType, on_delete=models.CASCADE)
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)
    url = models.URLField(max_length=URL_LENGTH, null=True)
    name = models.CharField(max_length=NAME_LENGTH)

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
    internal_id = models.CharField(max_length=150, null=True, help_text="The ID of the compound for internal use")
    operator = models.CharField(max_length=5, default="NA", help_text="An operator, like > < or =")

    def __str__(self) -> str:
        return f"{self.source}"

    def __repr__(self) -> str:
        return "<ActivityPoint %r %r %r %r %r %r>" % (self.id, self.source, self.target_id, self.activity, self.cmpd_id, self.units)

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
    description = models.CharField(max_length=255, default='',
                                   help_text='A short user-defined description for the project')
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, null=True, on_delete=models.CASCADE,
                                help_text='Foreign Key link to the relevant project'
                                          ' (optional for legacy reasons)')
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    tags = models.TextField(default='[]',
                            help_text='A comma-separated list of user-defined tags'
                                      ' used for searching and tagging projects')

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<SessionProject %r %r %r %r>" % (self.id, self.title, self.target, self.project)

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

    def __str__(self) -> str:
        return f"{self.author}"

    def __repr__(self) -> str:
        return "<SessionActions %r %r %r>" % (self.id, self.author, self.session_project)

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
        (MANUAL, 'MANUAL')  # Manual generated by user action
    )
    id = models.AutoField(primary_key=True)
    type = models.CharField(choices=SNAPSHOT_TYPE, default=INIT, max_length=8)
    title = models.CharField(max_length=255)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    description = models.CharField(max_length=255, default='')
    created = models.DateTimeField(default=timezone.now)
    data = models.TextField(help_text='JSON data that is passed from the front-end'
                                      ' describing what to load into the react components'
                                      ' to reproduce the session state')
    session_project = models.ForeignKey(SessionProject, null=True, on_delete=models.CASCADE)
    parent = models.ForeignKey('self', models.DO_NOTHING, blank=True, null=True, related_name='children',
                               help_text='Another Snapshot instance describing the current Snapshot parent')
    additional_info = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                       help_text='Optional JSON field containing name/value pairs for future use')

    def __str__(self) -> str:
        return f"{self.title}"

    def __repr__(self) -> str:
        return "<Snapshot %r %r %r %r>" % (self.id, self.title, self.type, self.author)

    class Meta:
        managed = True
        db_table = 'viewer_snapshot'


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
    session_project = models.ForeignKey(SessionProject, null=True, on_delete=models.CASCADE)
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE)
    last_update_date = models.DateTimeField(default=timezone.now)
    actions = models.JSONField(encoder=DjangoJSONEncoder)

    def __str__(self) -> str:
        return f"{self.author}"

    def __repr__(self) -> str:
        return "<SnapshotActions %r %r %r %r>" % (self.id, self.author, self.session_project, self.snapshot)

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
    compounds = models.ManyToManyField(Compound,
                                       help_text="The compounds that are in the design set")
    set_name = models.CharField(max_length=50)
    set_type = models.CharField(max_length=100, choices=SET_TYPE, default=USR)
    set_description = models.TextField(max_length=1000, blank=True, null=True)

    def __str__(self) -> str:
        return f"{self.set_name}"

    def __repr__(self) -> str:
        return "<DesignSet %r %r %r>" % (self.id, self.set_name, self.set_type)


class ComputedSetSubmitter(models.Model):
    name = models.CharField(max_length=50)
    email = models.CharField(max_length=100)
    institution = models.CharField(max_length=50,
                                   help_text="The institution or organizational affiliation"
                                             " of the compound set submitter")
    generation_date = models.DateField()
    method = models.CharField(max_length=50,
                              help_text="A name for the method that was used"
                                        " to produce the uploaded data")

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<ComputedSetSubmitter %r %r %r>" % (self.id, self.name, self.email)

    class Meta:
        unique_together = (("name", "method"),)


class CSetKeys(models.Model):
    """Used for authentication when uploading computed sets -
    each user is given an upload key associated with their email address in the form
    of a uuid. This is entered on the computed set upload page to allow a user upload.
    """
    user = models.CharField(max_length=50, default='User', editable=False)
    uuid = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True,
                            help_text="Unique key for the user")

    def __str__(self) -> str:
        return f"{self.user}"

    def __repr__(self) -> str:
        return "<CSetKeys %r %r %r>" % (self.id, self.user, self.uuid)


# computed sets = sets of poses calculated computationally
class ComputedSet(models.Model):
    """Computed sets - sets of 3D poses of molecules calculated computationally
    and uploaded by a user
    """
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
        (REVOKED, 'REVOKED')
    )

    name = models.CharField(max_length=50, unique=True, primary_key=True)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE)
    submitted_sdf = models.FileField(upload_to="compound_sets/", max_length=255,
                                     help_text="The sdf file containing the computed set")
    spec_version = models.FloatField(help_text="The version of the sdf file format specification")
    method_url = models.TextField(max_length=1000, null=True,
                                  help_text="A url linking to a write-up of the methodology used to create the"
                                            " computed set")
    submitter = models.ForeignKey(ComputedSetSubmitter, null=True, on_delete=models.CASCADE)
    unique_name = models.CharField(max_length=101,
                                   help_text="A unique name for the computed set,"
                                             " generated but the custom save() method")

    owner_user = models.ForeignKey(User, on_delete=models.CASCADE,
                                   default=settings.ANONYMOUS_USER)

    # The following fields will be used to track the computed set upload
    upload_task_id = models.CharField(null=True, max_length=50,
                                      help_text="The task ID of the upload Celery task")
    upload_status = models.CharField(choices=STATUS, null=True, max_length=7,
                                     help_text="Status of the upload. Only be updated at the end of the process")
    upload_progress = models.DecimalField(null=True, max_digits=5, decimal_places=2,
                                          help_text="Intended to be used as an indication of progress (0 to 100%)")
    upload_datetime = models.DateTimeField(null=True, help_text="The datetime the upload was completed")

    def save(self, **kwargs):
        """Custom save method for the ComputedSet model,
        including the method to auto-generate the unique name:
        unique_name = "".join(self.submitter.name.split()) + '-' + "".join(self.submitter.method.split())
        """
        if not self.submitter:
            super(ComputedSet, self).save()
        if not self.unique_name:
            unique_name = "".join(self.submitter.name.split()) + '-' + "".join(self.submitter.method.split())
            self.unique_name = unique_name
        super(ComputedSet, self).save()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<ComputedSet %r %r %r>" % (self.id, self.name, self.target)


class ComputedMolecule(models.Model):
    """The 3D information for a computed set molecule
    """
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)
    sdf_info = models.TextField(help_text="The 3D coordinates for the molecule")
    computed_set = models.ForeignKey(ComputedSet, on_delete=models.CASCADE)
    name = models.CharField(max_length=50)
    smiles = models.CharField(max_length=255)
    pdb = models.ForeignKey(Protein, on_delete=models.PROTECT, null=True)
    computed_inspirations = models.ManyToManyField(Molecule, blank=True,
                                                   help_text="Existing fragalysis molecules that were inspirations"
                                                             " in the design/calculation of the molecule")

    def __str__(self) -> str:
        return f"{self.smiles}"

    def __repr__(self) -> str:
        return "<ComputedMolecule %r %r %r %r>" % (self.id, self.smiles, self.name, self.compound)



class ScoreDescription(models.Model):
    """The names and descriptions of scores that the user uploads with each computed set molecule.
    """
    computed_set = models.ForeignKey(ComputedSet, null=True, on_delete=models.CASCADE)
    name = models.CharField(max_length=50, help_text="A name for this score")
    description = models.TextField(help_text="A description of this score,"
                                             " which should describe how to interpret it")

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<ScoreDescription %r %r>" % (self.id, self.name)


class NumericalScoreValues(models.Model):
    """The values of numerical scores that the user uploads with each computed set molecule.
    """
    score = models.ForeignKey(ScoreDescription, on_delete=models.CASCADE)
    value = models.FloatField()
    compound = models.ForeignKey(ComputedMolecule, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.score}"

    def __repr__(self) -> str:
        return "<NumericalScoreValues %r %r %r %r>" % (self.id, self.score, self.value, self.compound)



class TextScoreValues(models.Model):
    """The values of text scores that the user uploads with each computed set molecule.
    """
    score = models.ForeignKey(ScoreDescription, on_delete=models.CASCADE)
    value = models.TextField(max_length=500)
    compound = models.ForeignKey(ComputedMolecule, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.score}"

    def __repr__(self) -> str:
        return "<TextScoreValues %r %r %r %r>" % (self.id, self.score, self.value, self.compound)


class File(models.Model):
    file = models.FileField(blank=False)

    def __str__(self):
        return self.file.name

    def __repr__(self) -> str:
        return "<File %r %r>" % (self.id, self.file.name)


class DiscourseCategory(models.Model):
    """Discourse Subcategory references for Fragalysis - initially Targets
    """
    category_name = models.CharField(max_length=200, unique=True,
                                     help_text="The name of the (sub)category within Discourse. It must be unique")
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    discourse_category_id = models.IntegerField(help_text="The Discourse categoryID."
                                                          " Returned when the category was created")

    def __str__(self):
        return self.author

    def __repr__(self) -> str:
        return "<DiscourseCategory %r %r %r>" % (self.id, self.author, self.category_name)

    class Meta:
        db_table = 'viewer_discoursecategory'


class DiscourseTopic(models.Model):
    """Discourse Topic references for Fragalysis - initially Targets
    """
    topic_title = models.CharField(max_length=200, unique=True,
                                   validators=
                                   [MinLengthValidator(15,'Discourse Topic Title must be longer than 15 characters')],
                                   help_text="The title of the (sub)category within Discourse."
                                             " It must be unique within Discourse")
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    discourse_topic_id = models.IntegerField()

    def __str__(self):
        return self.author

    def __repr__(self) -> str:
        return "<DiscourseTopic %r %r %r>" % (self.id, self.author, self.topic_title)

    class Meta:
        db_table = 'viewer_discoursetopic'


class DownloadLinks(models.Model):
    """Searches made with the download_structures api.
    """
    file_url = models.CharField(max_length=200, unique=True, db_index=True,
                                help_text="Contains the complete link to the zip file including the uuid")
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE, db_index=True)
    proteins = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                help_text="Contains a sorted list of the protein codes in the search")
    protein_params = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                      help_text="Contains a sorted list of parameters used to create the zip file")
    other_params = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                    help_text="Contains a sorted list of parameters used to create the zip file")
    static_link = models.BooleanField(default=False,
                                      help_text="This preserves the proteins from the previous search")
    zip_contents = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                    help_text="For static files, this field contains the contents of the zip so that"
                                              " it can be reconstructed with the same file-links that it had previously."
                                              " For dynamic files, the zip is reconstructed from the search")
    create_date = models.DateTimeField()
    keep_zip_until = models.DateTimeField(db_index=True,
                                          help_text="The datetime when the tag was created"
                                                    " plus the retention time"
                                                    " (1 hour at the time of writing)")
    zip_file = models.BooleanField(default=False)
    original_search = models.JSONField(encoder=DjangoJSONEncoder, null=True)

    def __str__(self):
        return self.file_url

    def __repr__(self) -> str:
        return "<DownloadLinks %r %r %r %r>" % (self.id, self.file_url, self.user, self.target)

    class Meta:
        db_table = 'viewer_downloadlinks'


class TagCategory(models.Model):
    category = models.CharField(max_length=50, unique=True, help_text="The name of the tag category")
    colour = models.CharField(max_length=20, null=True, help_text="Expected to be an RGB string")
    description = models.CharField(max_length=200, null=True)

    def __str__(self):
        return self.category

    def __repr__(self) -> str:
        return "<TagCategory %r %r>" % (self.id, self.category)

    class Meta:
        db_table = 'viewer_tagcategory'


class Tag(models.Model):
    tag = models.CharField(max_length=200, help_text="The (unique) name of the tag")
    category = models.ForeignKey(TagCategory, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    create_date = models.DateTimeField(default=timezone.now)
    colour = models.CharField(max_length=20, null=True, help_text="Expected to be an RGB string")
    discourse_url = models.TextField(max_length=1000, null=True)
    help_text = models.TextField(null=True)
    additional_info = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                       help_text="Optional JSON field containing name/value pairs for future use")

    def __str__(self) -> str:
        return f"{self.tag}"

    def __repr__(self) -> str:
        return "<Tag %r %r %r %r %r>" % (self.id, self.tag, self.category, self.target, self.user)

    class Meta:
        abstract = True
        unique_together = ('tag', 'target',)


class MoleculeTag(Tag):
    molecules = models.ManyToManyField(Molecule, blank=True)
    mol_group = models.ForeignKey("scoring.MolGroup", null=True, blank=True,
                                  on_delete=models.SET_NULL)

    def __str__(self) -> str:
        return f"{self.mol_group}"

    def __repr__(self) -> str:
        return "<MoleculeTag %r %r>" % (self.id, self.mol_group)


class SessionProjectTag(Tag):
    """Data for SessionProjectTag(s) inherited from Tag.
    """
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
        (REVOKED, 'REVOKED')
    )
    id = models.AutoField(primary_key=True)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE, db_index=True)
    squonk_project = models.CharField(max_length=200, null=True)
    proteins = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    # Not used in phase 1
    compounds = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    transfer_spec = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    transfer_task_id = models.CharField(null=True, max_length=50)
    transfer_status = models.CharField(choices=STATUS, default=PENDING, max_length=7)
    transfer_progress = models.DecimalField(null=True, max_digits=5, decimal_places=2,
                                            help_text="Intended to be used as an indication of progress (0 to 100%)")
    transfer_datetime = models.DateTimeField(null=True, help_text="The datetime the transfer was completed")

    def __str__(self) -> str:
        return f"{self.user}"

    def __repr__(self) -> str:
        return "<JobFileTransfer %r %r %r %r %r>" % (self.id, self.user, self.snapshot, self.target, self.squonk_project)

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
        (SUCCESS, 'SUCCESS'),  # Job completed successfully in squonk (updated by squonk)
        (FAILURE, 'FAILURE'),  # Job failed in squonk (updated by squonk)
        (RETRY, 'RETRY'),  # Job status in squonk (updated by squonk)
        (REVOKED, 'REVOKED')  # Job status in squonk (updated by squonk)
    )
    UPLOAD_STATUS = (
        (PENDING, 'PENDING'),  # Initial state when upload queued
        (STARTED, 'STARTED'),  # Upload job started
        (SUCCESS, 'SUCCESS'),  # Upload job successful
        (FAILURE, 'FAILURE'),  # Upload job failed
        (RETRY, 'RETRY'),
        (REVOKED, 'REVOKED')
    )
    id = models.AutoField(primary_key=True)
    squonk_job_name = models.CharField(max_length=200, null=True)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE,
                                 help_text="The snapshot the file transfer is part of")
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE,
                               db_index=True)
    project = models.ForeignKey(Project, null=True, on_delete=models.CASCADE)
    squonk_project = models.CharField(max_length=200, null=True,
                                      help_text="The name of a project that has been created in Squonk"
                                                " that the files will be transferred to")
    squonk_job_spec = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    job_start_datetime = models.DateTimeField(null=True)
    job_finish_datetime = models.DateTimeField(null=True,
                                               help_text="The datetime when the Squonk Job has finished,"
                                                         " populated by information in the Squonk callback."
                                                         " If this is not set you can assume the JOb is still"
                                                         " running. When it is set the job_status filed will be"
                                                         " updated (to SUCCESS or FAILURE). If automatic upload"
                                                         " follows an upload_task_id wil be set and you can"
                                                         " monitor upload_status for a status of the upload")
    job_status = models.CharField(choices=SQUONK_STATUS, default=PENDING, max_length=7,
                                  help_text="The status of the Squonk job, e.g. 'PENDING'"
                                            " Will be modified by Squonk through the callback URL")
    job_status_datetime = models.DateTimeField(null=True)
    squonk_job_info = models.JSONField(encoder=DjangoJSONEncoder, null=True,
                                       help_text="squonk_job_info is a copy of the response from DmApi.start_job_instance()."
                                                 " It's an instance of a DmApiRv object (a namedtuple)"
                                                 " that contains a 'success' (boolean) and 'msg' (the DmApi response's resp.json())."
                                                 " For us this will contain a 'task_id', 'instance_id' and 'callback_token'."
                                                 " The content will be a list with index '0' that's the value of the DmApiRv"
                                                 " 'success' variable and, at index '1', the original response message json()."
                                                 " The Job callback token will be squonk_job_info[1]['callback_token']")
    squonk_url_ext = models.CharField(max_length=200, null=True,
                                      help_text="a Squonk UI URL to obtain information about the running instance."
                                                " It's essentially the Squonk URL with the instance ID appended.")
    code = models.UUIDField(default=uuid.uuid4, editable=False, unique=True,
                            help_text="A UUID generated by Fragalysis and passed to Squonk as part of a callback URL."
                                      " This value is used to uniquely identify the HJob in Squonk and is passed back"
                                      " by squonk to provide context in calls to the JobCallBackView")
    upload_task_id = models.CharField(null=True, max_length=50,
                                      help_text="Celery task ID for results upload task."
                                                " Set when the Job completes and an automated upload follows")
    upload_status = models.CharField(choices=UPLOAD_STATUS, default=PENDING, max_length=7,
                                     null=True, help_text="Status of upload task")
    computed_set = models.ForeignKey(ComputedSet, on_delete=models.CASCADE, null=True)

    def __str__(self) -> str:
        return f"{self.user}"

    def __repr__(self) -> str:
        return "<JobRequest %r %r %r %r %r %r>" % (self.id, self.user, self.squonk_job_name, self.snapshot, self.target, self.squonk_project)

    class Meta:
        db_table = 'viewer_jobrequest'


class JobOverride(models.Model):
    override = models.JSONField(encoder=DjangoJSONEncoder)
    author = models.ForeignKey(User, null=True, on_delete=models.SET_NULL,
                               help_text="The user that uploaded the override")

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
    uuid = models.TextField(max_length=40,
                            help_text="A Squonk2 Account Server (AS) Organisation UUID."
                                      " A fixed length string consisting of 'org-' followed by a uuid4 value,"
                                      " e.g. 'org-54260047-183b-42e8-9658-385a1e1bd236'")
    name = models.TextField(max_length=80,
                            help_text="The name of the Squonk2 Organisation UUID (obtained form the AS)")
    as_url = models.URLField()
    as_version = models.TextField()

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Squonk2Org %r %r %r>" % (self.id, self.name, self.uuid)


class Squonk2Unit(models.Model):
    """Squonk2 Unit (UUIDs). Managed by the Squonk2Agent class.
    """
    uuid = models.TextField(max_length=41,
                            help_text="A Squonk2 Account Server (AS) Unit UUID."
                                      " A fixed length string consisting of 'unit-' followed by a uuid4 value,"
                                      " e.g. 'unit-54260047-183b-42e8-9658-385a1e1bd236'")
    name = models.TextField(help_text="The name used to create the Squonk2 Unit UUID"
                                      " This is not limited by the actual name length imposed by the DM")
    organisation = models.ForeignKey(Squonk2Org, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Squonk2Unit %r %r %r>" % (self.id, self.name, self.uuid)


class Squonk2Project(models.Model):
    """Squonk2 Project and Product (UUIDs). Managed by the Squonk2Agent class.
    """
    uuid = models.TextField(max_length=44)
    name = models.TextField(help_text="The name of the Squonk2 Project UUID (obtained form the AS)")
    product_uuid = models.TextField(max_length=44,
                                    help_text="A Squonk2 Account Server (AS) Product UUID."
                                              " A fixed length string consisting of 'product-' followed by a uuid4 value,"
                                              " e.g. 'product-54260047-183b-42e8-9658-385a1e1bd236'")
    unit = models.ForeignKey(Squonk2Unit, on_delete=models.CASCADE)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<Squonk2Project %r %r %r %r %r>" % (self.id, self.name, self.uuid, self.product_uuid, self.unit)



class CanonSite(models.Model):
    name = models.TextField()
    quat_assembly = models.ForeignKey(QuatAssembly, null=True, on_delete=models.CASCADE)
    residues = models.JSONField(encoder=DjangoJSONEncoder)
    # TODO: missing in db, check if correct, (might be correct, but might not)
    ref_conf_site = models.OneToOneField("CanonSiteConf", null=True, on_delete=models.CASCADE)
    canon_site_id = models.IntegerField(null=False, help_text="canon_site id from YAML")

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CanonSite %r %r>" % (self.id, self.name)


class XtalformSite(models.Model):
    xtalform = models.ForeignKey(Xtalform, on_delete=models.CASCADE)
    canon_site = models.ForeignKey(CanonSite, on_delete=models.CASCADE)
    lig_chain = models.CharField(max_length=1)
    residues = models.JSONField(encoder=DjangoJSONEncoder)
    xtalform_site_id = models.IntegerField(null=False, help_text="xtalform site id from YAML")

    def __str__(self) -> str:
        return f"{self.xtalform_site_id}"

    def __repr__(self) -> str:
        return "<CanonSiteConf %r %r %r>" % (self.id, self.xtalform_site_id, self.xtalform)


class CanonSiteConf(models.Model):
    canon_site = models.ForeignKey(CanonSite, on_delete=models.CASCADE)
    # TODO: name not present in metadata atm
    name = models.TextField(null=True)
    ref_site_observation = models.OneToOneField("SiteObservation", null=True, on_delete=models.CASCADE)
    residues = models.JSONField(encoder=DjangoJSONEncoder)

    def __str__(self) -> str:
        return f"{self.name}"

    def __repr__(self) -> str:
        return "<CanonSiteConf %r %r %r>" % (self.id, self.name, self.canon_site)


class SiteObservation(models.Model):
    code = models.TextField()
    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)
    cmpd = models.ForeignKey(Compound, on_delete=models.CASCADE)
    xtalform_site = models.ForeignKey(XtalformSite, on_delete=models.CASCADE)
    canon_site_conf = models.ForeignKey(CanonSiteConf, on_delete=models.CASCADE)
    bound_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    apo_solv_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    apo_desolv_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    apo_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    xmap_2fofc_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    event_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    artefacts_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    pdb_header_file = models.FileField(upload_to="target_loader_data/", null=True, max_length=255)
    smiles = models.TextField()
    seq_id = models.IntegerField()
    chain_id = models.CharField(max_length=1)
    ligand_mol_file = models.TextField(null=True)

    history = HistoricalRecords()

    def __str__(self) -> str:
        return f"{self.code}"

    def __repr__(self) -> str:
        return "<SiteObservation %r %r %r %r>" % (self.id, self.code, self.experiment, self.compound)
