import uuid
from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone
from django.core.serializers.json import DjangoJSONEncoder
from django.core.validators import MinLengthValidator
from django.conf import settings

from shortuuid.django_fields import ShortUUIDField

from simple_history.models import HistoricalRecords

from viewer.target_set_config import get_mol_choices, get_prot_choices

class Project(models.Model):
    """Django model for holding information about a project. This is used on the Targets level, adding a new project for
    each target, and saving a list of users that can access a target during the authentication step

    Parameters
    ----------
    title: CharField
        The title of the project
    init_date: DateTimeField
        The date the project was initiated (autofield)
    user_id: ManyToManyField
        Links to the User model
    open_to_public: BooleanField
        True if open to the Public
    """
    # The title of the project_id -> userdefined
    title = models.CharField(max_length=200, unique=True)
    # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    # The users it's related to
    user_id = models.ManyToManyField(User)
    open_to_public = models.BooleanField(default=False)


class Target(models.Model):
    """Django model to define a Target - a protein.

    Parameters
    ----------
    title: CharField
        The name of the target
    init_date: DateTimeField
        The date the target was initiated (autofield)
    project_id: ManyToManyField
        Links targets to projects for authentication
    uniprot_id: Charfield
        Optional field where a uniprot id can be stored
    metadata: FileField
        Optional file upload defining metadata about the target - can be used to add custom site labels
    zip_archive: FileField
        Link to zip file created from targets uploaded with the loader
    default_squonk_project = CharField
        Contains the default Squonk project name that jobs will be run in
    upload_task_id: CharField
        Task id of upload celery task Note that if a resynchronisation is
        required this will be re-used.
    upload_status: CharField
        Identifies the status of the upload (note will only be updated at the end of the process.
    upload_progress: DecimalField
        Intended to be used as an indication of progress (0 to 100%)
    upload_datetime: DateTimeField
        The datetime the upload was completed.
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

    # The title of the project_id -> userdefined
    title = models.CharField(unique=True, max_length=200)
    # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    # A field to link projects and targets together
    project_id = models.ManyToManyField(Project)
    # Indicates the uniprot_id id for the target. Is a unique key
    uniprot_id = models.CharField(max_length=100, null=True)
    metadata = models.FileField(upload_to="metadata/", null=True, max_length=255)
    zip_archive = models.FileField(upload_to="archive/", null=True, max_length=255)
    default_squonk_project = models.CharField(max_length=200, null=True)
    # The following fields will be used to track the target upload
    upload_task_id = models.CharField(null=True, max_length=50)
    upload_status = models.CharField(choices=STATUS, null=True, max_length=7)
    upload_progress = models.DecimalField(null=True, max_digits=5, decimal_places=2)
    upload_datetime = models.DateTimeField(null=True)


class Protein(models.Model):
    """Django model for holding information about a protein. A protein is a protein structure which has a unique set of
    3D coordinates, rather than a target, which is a set of protein objects of the same protein. A Molecule object is
    also linked to a protein, so that a complete structure is comprised of the molecule and protein in separate parts in
    fragalysis.

    Parameters
    ----------
    code: CharField
        A unique name for the protein (e.g. NUDT5A-x0001_1)
    target_id: ForeignKey
        Foreign key linking the protein to it's target
    apo_holo: NullBooleanField
        0 for apo (ligand removed), 1 for holo (ligand in tact)
    prot_type: CharField
        protein type - from a pre-defined list and determined by file extension on upload (defined in
        loader.config.get_prot_choices):
            prot_choices = (
            (APO, "Apo", "_apo.pdb", "APO"),
            (STRIPPED, "Stripped", "_no_buffer_altlocs.pdb", "STRIPPED"),
            (TLEAPED, "Tleaped", "_tleap.pdb", "TLEAP"),
            (CHUNKED, "Chunked", "_chunk.pdb", "CHUNK"),
            (BOUND, "Bound", '_bound.pdb', "BOUND")
            )
    pdb_info: FileField
        File link to apo pdb structure - pdb file with ligand removed
    bound_info: FileField
        File link to bound state structure - same as apo pdb but with ligand in-tact
    cif_info: FileField
        File link to cif file for ligand restraints (optional)
    mtz_info: FileField
        File link to uploaded mtz file (optional)
    map_info: FileField
        File link to uploaded map file (optional)
    trans_matrix_info: FileField
        File link to uploaded transformation matrix file (optional)
    pdb_header_info: FileField
        File link to uploaded _header.pdb file (optional)
    apo_desolve_info: FileField
        File link to uploaded _apo-desolv.pdb file (optional)
    aligned: NullBooleanField
        Bool - 1 if aligned, 0 if not
    aligned_to: ForeignKey (self)
        Foreign key to another instance of a protein to which this protein is aligned (optional)
    has_eds: NullBooleanField
        Bool - 1 if has ED, 0 it not (optional)
    """
    # code for this protein (e.g. NUDT5A-x0001_1A)
    code = models.CharField(max_length=50, db_index=True)
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    apo_holo = models.NullBooleanField()
    # Set the groups types
    prot_choices, default_prot = get_prot_choices()
    prot_type = models.CharField(
        choices=prot_choices, default=default_prot, max_length=2
    )
    pdb_info = models.FileField(upload_to="pdbs/", null=True, max_length=255)
    bound_info = models.FileField(upload_to="bound/", null=True, max_length=255)
    cif_info = models.FileField(upload_to="cifs/", null=True, max_length=255)
    mtz_info = models.FileField(upload_to="mtzs/", null=True, max_length=255)
    map_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    sigmaa_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    diff_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    event_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    trans_matrix_info = models.FileField(upload_to="trans/", null=True, max_length=255)
    pdb_header_info = models.FileField(upload_to="pdbs/", null=True, max_length=255)
    apo_desolve_info = models.FileField(upload_to="pdbs/", null=True, max_length=255)
    aligned = models.NullBooleanField()
    aligned_to = models.ForeignKey("self", null=True, on_delete=models.CASCADE)
    has_eds = models.NullBooleanField()

    class Meta:
        unique_together = ("code", "target_id", "prot_type")


class Compound(models.Model):
    """Django model for holding information about a compound, which is a unique 2D molecule

    Parameters
    ----------
    inchi: CharField
        The inchi key of the compound
    long_inchi: TextField
        For historical reasons, the inchifield cannot be removed, but has a max limit of 255 characters. If the inchi
        key for a compound is longer than this, it is stored in the long_inchi field, and the inchi key is concatenated
        to the first 255 characters. (optional)
        TODO: Use a method to get the inchi key
    smiles: Charfield
        The SMILES string representation of the compound
    current_identifier: Charfield
        The identifier for this compound that is used in Fragalysis to represent it's 3D molecule (optional)
    all_identifiers: TextField
        A comma separated list of all identifiers that have been used in the past to represent this 2D compound
    project_id: ManyToManyField
        Many to Many foreign key relationship to all projects that this compound is associated to (prevents duplication
        of compounds across multiple targets)
    mol_log_p: FloatField
        Computed LogP value (from rdkit)
    mol_wt: FloatField
        Computed molecular weight (Da) (from rdkit)
    tpsa: FloatField
        Computed Topological Polar Surface Area (from rdkit)
    heavy_atom_count: IntegerField
        Computed heavy (non-hydrogen) atom count (from rdkit)
    heavy_atom_mol_wt: FloatField
        Computed molecular weight of all heavy (non-hydrogen) atoms (from rdkit)
    nhoh_count: IntegerField
        Computed number of hydroxylamine groups (from rdkit)
    no_count: IntegerField
        Computed number of nitroso groups (from rdkit)
    num_h_acceptors: IntegerField
        Computed number of hydrogen-bond acceptor groups (from rdkit)
    num_h_donors: IntegerField
        Computed number of hydrogen-bond donor groups (from rdkit)
    num_het_atoms: IntegerField
        Computed number of heterogeneous atoms (from rdkit)
    num_rot_bonds: IntegerField
        Computed number of rotatable bonds (from rdkit)
    num_val_electrons: IntegerField
        Computed number of valence electrons (from rdkit)
    ring_count: IntegerField
        Computed number of rings in the molecule (from rdkit)
    inspirations: ManyToManyField
        Foreign key link to any number of 3D Molecules that inspired the design of this compound
    description: TextField
        A description of the compound added by a user (optional)
    comments: TextField
        A free-text comments field (optional)
    """

    # Character attributes
    inchi = models.CharField(max_length=255, unique=False, db_index=True)
    long_inchi = models.TextField(max_length=1000, blank=True, null=True)
    smiles = models.CharField(max_length=255, db_index=True)
    current_identifier = models.CharField(max_length=255, db_index=True, blank=True, null=True)
    all_identifiers = models.TextField(blank=True, null=True)
    # A link to the related project
    project_id = models.ManyToManyField(Project)
    # Float attributes
    mol_log_p = models.FloatField()
    mol_wt = models.FloatField()
    tpsa = models.FloatField()
    # Integer attributes
    heavy_atom_count = models.IntegerField()
    heavy_atom_mol_wt = models.FloatField()
    nhoh_count = models.IntegerField()
    no_count = models.IntegerField()
    num_h_acceptors = models.IntegerField()
    num_h_donors = models.IntegerField()
    num_het_atoms = models.IntegerField()
    num_rot_bonds = models.IntegerField()
    num_val_electrons = models.IntegerField()
    ring_count = models.IntegerField()
    inspirations = models.ManyToManyField("Molecule", blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    comments = models.TextField(blank=True, null=True)

    class Meta:
        unique_together = ('inchi', 'long_inchi')


class Molecule(models.Model):
    """Django model for holding information about a Molecule. A molecule is linked to a compound, and represents the
    unique 3D coordinates for that molecule. Note a compound can be linked to multiple molecules, each with a different
    3D structure, but the same 2D structure

    Parameters
    ----------
    smiles: CharField
        The smiles string of the molecule. (optional) TODO: Check if this is needed...
    lig_id: CharField
        The ligand ID from the pdb file it originated from (e.g. LIG 1 A) (optional)
    chain_id: CharField
        The chain from the protein structure pdb file that this ligand is from (e.g. A) (optional)
    mol_type: CharField
        molecule type - from a pre-defined list and determined by file extension on upload (defined in
        loader.config.get_mol_choices):
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
    sdf_info: TextField
        The 3D coordinates for the molecule in MDL (mol file) format. Taken directly from the uploaded file
    rscc: FloatField
        The RSCC score of the molecule 3D coordinates vs. PANDDA event map (optional)
    occupancy: FloatField
        The occupancy (electron density) of the molecule (optional)
    x_com: FloatField
        x-coordinate for centre of mass (optional)
    y_com: FloatField
        y-coordinate for centre of mass (optional)
    z_com: FloatField
        z-coordinate for centre of mass (optional)
    rmsd: FloatField
        RMSD of this molecule compared to ? (optional - unused)
    prot_id: ForeignKey
        Foreign key link to the associated protein (apo) that this ligand was pulled from
    cmpd_id: ForeignKey
        Foreign key link to the associated 2D compound
    sdf_file: FileField
        File link to uploaded sdf file (optional)
    history: HistoricalRecords
        Tracks the changes made to an instance of this model over time

    """
    # Character attributes
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    lig_id = models.CharField(max_length=5, null=True)
    chain_id = models.CharField(max_length=1, null=True)
    # The type of map
    # Set the groups types
    mol_choices, default_mol = get_mol_choices()
    mol_type = models.CharField(choices=mol_choices, default=default_mol, max_length=2)
    # Textfield
    sdf_info = models.TextField(null=True)
    # Float attributes
    rscc = models.FloatField(null=True)
    occupancy = models.FloatField(null=True)
    x_com = models.FloatField(null=True)
    y_com = models.FloatField(null=True)
    z_com = models.FloatField(null=True)
    rmsd = models.FloatField(null=True)
    # Foreign key relations
    prot_id = models.ForeignKey(Protein, on_delete=models.CASCADE)
    cmpd_id = models.ForeignKey(Compound, on_delete=models.CASCADE)
    sdf_file = models.FileField(upload_to="sdfs/", null=True, max_length=255)
    history = HistoricalRecords()

    # Unique constraints
    class Meta:
        unique_together = ("prot_id", "cmpd_id", "mol_type")


class ActivityPoint(models.Model):
    """Django model for holding information about the activity of a compound - currently unused

    Parameters
    ----------
    source: CharField
        The source of the activity data
    target_id: ForeignKey
        Foreign key link to the relevant target
    cmpd_id: ForeignKey
        ForeignKey link to the relevant compound
    activity: FloatField
        Measured -log(10) activity
    units: Charfield
        The units (e.g. uM or whatever) TODO: Convert to choices
    confidence: IntegerField
        The confidence level set for the activity measurement (optional)
    internal_id: CharField
        The ID of the compound for internal use
    operator: Charfield
        The operator > < or = (= is default)
    """
    # This should encompass the activity type too
    source = models.CharField(max_length=50, null=True, db_index=True)
    # Foreign key to Target object
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    # Foreign key to compound object
    cmpd_id = models.ForeignKey(Compound, on_delete=models.CASCADE)
    # Measured -log(10) activity
    activity = models.FloatField(db_index=True)
    # The units (uM or whatever). Should be more option
    units = models.CharField(max_length=50)
    # The confidence score
    confidence = models.IntegerField(null=True, db_index=True)
    # The ID of the compound for internal use
    internal_id = models.CharField(max_length=150, null=True)
    # The operator > < or = (= is default)
    operator = models.CharField(max_length=5, default="NA")

    class Meta:
        unique_together = ("target_id", "activity", "cmpd_id", "units")

# Start of Action Types
class ActionType(models.Model):
    """Django model for holding types of the actions that may be performed by a user during a session or snapshot.

    Parameters
    ----------
    id: Autofield
        Auto-created id for the action type
    description: CharField
        The description of the action type
    active: NullBooleanField
        Bool - True if action type is active, False if not (historical, not active yet, not important)
    activation_date: DateTimeField
        The datetime the action type became active (set with a change of status)
    """
    id = models.AutoField(primary_key=True)
    description = models.CharField(max_length=200, default='')
    active = models.BooleanField(default=False)
    activation_date = models.DateTimeField(default=timezone.now)

    class Meta:
        db_table = 'viewer_actiontype'


# Start of Session Project
class SessionProject(models.Model):
    """Django model for holding information about a fragalysis user Session Project
    - a set of sessions saved by a user that belong to a Target and Project.

    Parameters
    ----------
    title: CharField
        The title of the project
    init_date: DateTimeField
        The date the project was initiated (autofield)
    description: Charfield
        A short user-defined description for the project
    target: ForeignKey
        Foreign Key link to the relevant project target
    project: ForeignKey
        Foreign Key link to the relevant project (optional for legacy reasons)
    author: ForeignKey
        A link to the user that created the project
    tags: TextField
        A comma separated list of user-defined tags - for searching and tagging projects
    """
    title = models.CharField(max_length=200)
    init_date = models.DateTimeField(default=timezone.now)
    description = models.CharField(max_length=255, default='')
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    project = models.ForeignKey(Project, null=True, on_delete=models.CASCADE)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    tags = models.TextField(default='[]')

    class Meta:
        db_table = 'viewer_sessionproject'

class SessionActions(models.Model):
    """Django model for holding the user actions related to a particular session_project.

    Parameters
    ----------
    id: Autofield
        Auto-created id for the session actions table, initially should be a 1 to 1 with session project.
    author_id: ForeignKey
        Foreign key link to the id of the user that created the session_project
    session_project: ForeignKey
        A foreign key link to the relevant session_project (required)
    last_update_date: DateTimeField
        Timestamp for when the action list was generated or updated
    actions : JSONField
        A JSON field containing types of actions related to the session_project.
        The list elements are (at the time of writing) :
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

    class Meta:
        db_table = 'viewer_sessionactions'


class Snapshot(models.Model):
    """Django model for storing information describing a static snapshot of a fragalysis page - multiple snapshots make
    up a project

    Parameters
    ----------
    id: Autofield
        Auto-created id for the session, used in url accession
    type: CharField
        Describes the session type:
            SNAPSHOT_TYPE = (
            (INIT, "INIT"),  # Initial snapshot generated by system
            (AUTO, 'AUTO'),  # Automatic generated by system
            (MANUAL, 'MANUAL')  # Manual generated by user action
        )
    title: Charfield
        The title of the snapshot
    author: ForeignKey
        Foreign key link to the user that created the snapshot
    description: Charfield
        Short user-provided description for the snapshot
    created: DateTimeField
        Auto-created timestamp for when the snapshot was created
    data: TextField
        Field to hold the json data that is passed from the front-end describing what to load into the react components
        to reproduce the session state
    session_project: ForeignKey
        If the snapshot is part of a project, a foreign key link to the relevant project (optional)
    parent: ForeignKey(self)
        Foreign key link to another Snapshot instance describing the current Snapshot parent (optional)
    additional_info: JSONField
        Optional JSON field containing name/value pairs for future use
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
    data = models.TextField()
    session_project = models.ForeignKey(SessionProject, null=True, on_delete=models.CASCADE)
    parent = models.ForeignKey('self', models.DO_NOTHING, blank=True, null=True, related_name='children')
    additional_info = models.JSONField(encoder=DjangoJSONEncoder, null=True)

    class Meta:
        managed = True
        db_table = 'viewer_snapshot'


class SnapshotActions(models.Model):
    """Django model for holding the user actions leading to a particular snapshot or idea.

    Parameters
    ----------
    id: Autofield
        Auto-created id for the session actions table, initially should be a 1 to 1 with snapshot.
    author_id: ForeignKey
        Foreign key link to the id of the user that created the snapshot
    session_project: ForeignKey
        If the snapshot is part of a project, a foreign key link to the relevant project (optional)
    snapshot: ForeignKey
        A foreign key link to the relevant snapshot (required)
    last_update_date: DateTimeField
        Timestamp for when the action list was generated or updated
    actions : JSONField
        A JSON field containing types of actions made leading to the snapshot.
        The list elements are (at the time of writing) :
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

    class Meta:
        db_table = 'viewer_snapshotactions'


# End of Session Project


# Start of compound sets
# Design sets = 2D compounds that have been designed but do not yet have a 3D structure
class DesignSet(models.Model):
    """Django model for holding information about a design sets - sets of 2D compounds that have been designed but do
    not yet have a 3D structure - unused

    Parameters
    ----------
    compounds: ManyToManyField
        Foreign key links to all 2D compounds that are in the design set
    set_name: CharField
        A user-defined name for the design set
    set_type: Charfield
        The type of design set, from a pre-defined set of choices:
            SET_TYPE = (
                    (LIB, "library"),  # library - e.g. DSiPoised
                    (FUP, 'follow-up'),  # follow-up - e.g. purchased compounds
                    (USR, 'user-submitted'),  # user submitted - can be submitted by anyone
                    (ENM, 'enumerated'),  # enumerated - e.g. similarity search or something
                )
    set_description: TextField
        A user-input description of the design set (optional)
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
    compounds = models.ManyToManyField(Compound)
    set_name = models.CharField(max_length=50)
    set_type = models.CharField(max_length=100, choices=SET_TYPE, default=USR)
    set_description = models.TextField(max_length=1000, blank=True, null=True)


class ComputedSetSubmitter(models.Model):
    """Django model for holding information about the submitter of a computed set

    Parameters
    ----------
    name: CharField
        The name of the computed set submitter
    email: CharField
        The email address of the computed set submitter
    institution: Charfield
        The institution or organizational affiliation of the compound set submitter
    generation_date: DateField
        The date that the uploaded data was generated on
    method: Charfield
        A name for the method that was used to produce the uploaded data
    """

    name = models.CharField(max_length=50, null=False)
    email = models.CharField(max_length=100, null=False)
    institution = models.CharField(max_length=50, null=False)
    generation_date = models.DateField()
    method = models.CharField(max_length=50, null=False)

    class Meta:
        unique_together = (("name", "method"),)


class CSetKeys(models.Model):
    """Django model for authentication when uploading computed sets - each user is given an upload key associated with
    their email address in the form of a uuid. This is entered on the computed set upload page to allow a user upload.

    Parameters
    ----------
    user: CharField
        User email address
    uuid: UUIDField
        Auto-generated unique uuid4 string for the user
    """
    user = models.CharField(max_length=50, default='User', editable=False)
    uuid = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True)


# computed sets = sets of poses calculated computationally
class ComputedSet(models.Model):
    """Django model holding information about computed sets - sets of 3D poses of molecules calculated computationally
    and uploaded by a user

    Parameters
    ----------
    name: CharField
        A unique name for the computed set
    target: ForeignKey
        Foreign key link to the relevant target
    submitted_sdf: FileField
        File link to a stored version of the sdf file that the user uploaded
    spec_version: FloatField
        Version number of the sdf file format specification for upload of computed sets
    method_url: TextField
        A url linking to a write-up of the methodology used to create a computed set
    submitter: ForeignKey
        Foreign key link to the submitter information
    unique_name: CharField
        Auto-generated unique name for a computed set
    owner_user: ForeignKey
        A link to the user that created the Computed Set
    upload_task_id: CharField
        Task id of upload celery task Note that if a resynchronisation is
        required this will be re-used.
    upload_status: CharField
        Identifies the status of the upload (note will only be updated at the end of the process.
    upload_progress: DecimalField
        Intended to be used as an indication of progress (0 to 100%)
    upload_datetime: DateTimeField
        The datetime the upload was completed.
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

    # a (unique) name for this compound set
    name = models.CharField(max_length=50, unique=True, primary_key=True)
    # target that this compound set belongs to
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE)
    # link to the submitted sdf file
    submitted_sdf = models.FileField(upload_to="compound_sets/", null=False, max_length=255)
    # file format specification version
    spec_version = models.FloatField(null=False)
    method_url = models.TextField(max_length=1000, null=True)
    submitter = models.ForeignKey(ComputedSetSubmitter, null=True, on_delete=models.CASCADE)
    unique_name = models.CharField(max_length=101, null=False)

    owner_user = models.ForeignKey(User, null=False, on_delete=models.CASCADE,
                                      default=settings.ANONYMOUS_USER)

    # The following fields will be used to track the computed set upload
    upload_task_id = models.CharField(null=True, max_length=50)
    upload_status = models.CharField(choices=STATUS, null=True, max_length=7)
    upload_progress = models.DecimalField(null=True, max_digits=5, decimal_places=2)
    upload_datetime = models.DateTimeField(null=True)

    # Check if needed? Rachael still having a look.
    # design_set = models.ForeignKey(DesignSet, null=False, blank=False)

    def save(self, **kwargs):
        """Custom save method for the ComputedSet model, including the method to auto-generate the unique name:
            unique_name = "".join(self.submitter.name.split()) + '-' + "".join(self.submitter.method.split())
        """
        if not self.submitter:
            super(ComputedSet, self).save()
        if not self.unique_name:
            unique_name = "".join(self.submitter.name.split()) + '-' + "".join(self.submitter.method.split())
            self.unique_name = unique_name
        super(ComputedSet, self).save()


class ComputedMolecule(models.Model):
    """Django model to hold the 3D information for a computed set molecule

    Parameters
    ----------
    compound: ForeignKey
       Foreign key link to the molecule 2D information
    sdf_info: TextField
       The 3D coordinates for the molecule in MDL (mol file) format. Taken directly from the uploaded file
    computed_set: ForeignKey
        Foreign key link to the computed set that this molecule is a part of
    name: CharField
        A name for the molecule
    smiles: Charfield
        SMILES string for the molecule
    pdb_info: FileField
        A file link to a user-uploaded apo structure for this molecule, if an existing fragalysis protein was not used
        (optional)
    computed_inspirations: ManyToManyField
        Foreign key links to existing fragalysis molecules that were inspirations in the design/calculation of the
        molecule

    """
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)
    # the 3D coordinates of a computed molecule
    sdf_info = models.TextField(null=False)
    # a link to the compound set this molecule belongs to
    computed_set = models.ForeignKey(ComputedSet, on_delete=models.CASCADE)
    # a name for this compound
    name = models.CharField(max_length=50)
    # calculated smiles
    smiles = models.CharField(max_length=255)
    pdb = models.ForeignKey(Protein, on_delete=models.PROTECT, null=True)
    # the protein link
    # pdb_info = models.FileField(upload_to="pdbs/", null=False, max_length=255)
    # if we use our own method of calculating them
    computed_inspirations = models.ManyToManyField(Molecule, null=True, blank=True)

    # @property
    # def pdb_info(self):
    #     if self.pdb:
    #         return self.pdb.pdb_info
    #     else:
    #         return None


class ScoreDescription(models.Model):
    """Django model to store the names and descriptions of scores that the user uploads with each computed set molecule

    Parameters
    ----------
    computed_set: ForeignKey
        Foreign key link to the relevant computed set
    name: Charfield
        A name for the score
    description: TextField
        A description of the score, which should describe how to interpret it
    """
    # which compound set this score belongs to
    computed_set = models.ForeignKey(ComputedSet, null=True, on_delete=models.CASCADE)
    # a name for this score
    name = models.CharField(max_length=50)
    # a description for this score
    description = models.TextField(null=False)


class NumericalScoreValues(models.Model):
    """Django model to store the values of numerical scores that the user uploads with each computed set molecule

    Parameters
    ----------
    score: ForeignKey
        Foreign key link to the relevant score name and description
    value: FloatField
        The numerical value for the score
    compound: ForeignKey
        Foreign key link to the computed molecule that the score corresponds to
    """
    # a link to the score name and description
    score = models.ForeignKey(ScoreDescription, on_delete=models.CASCADE)
    # the specific value for this score
    value = models.FloatField(null=False)
    # the compound this score relates to
    compound = models.ForeignKey(ComputedMolecule, on_delete=models.CASCADE)


class TextScoreValues(models.Model):
    """Django model to store the values of text scores that the user uploads with each computed set molecule

       Parameters
       ----------
       score: ForeignKey
           Foreign key link to the relevant score name and description
       value: TextField
           The text value for the score
       compound: ForeignKey
           Foreign key link to the computed molecule that the score corresponds to
       """
    # a link to the score name and description
    score = models.ForeignKey(ScoreDescription, on_delete=models.CASCADE)
    # the specific value for this score
    value = models.TextField(max_length=500, null=False)
    # the compound this score relates to
    compound = models.ForeignKey(ComputedMolecule, on_delete=models.CASCADE)
# End of compound sets


class File(models.Model):
    file = models.FileField(blank=False, null=False)

    def __str__(self):
        return self.file.name


# Start of Discourse Translation Tables
class DiscourseCategory(models.Model):
    """Django model for holding Discourse Subcategory references for Fragalysis - initially Targets

    Parameters
    ----------
    category_name: CharField
        The name of the (sub)category within Discourse. It must be unique within Discourse
    author: ForeignKey
        A link to the user that created the category
    discourse_category_id: IntegerField
        The Discourse category_id returned when the category was created. Used when creating new topics.

    """
    category_name = models.CharField(max_length=200, unique=True)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    discourse_category_id = models.IntegerField()

    class Meta:
        db_table = 'viewer_discoursecategory'


class DiscourseTopic(models.Model):
    """Django model for holding Discourse Topic references for Fragalysis - initially Targets

    Parameters
    ----------
    topic_title: CharField
        The title of the sub)category within Discourse. It must be unique within Discourse.
    author: ForeignKey
        A link to the user that created the category
    discourse_topic_id: IntegerField
        The Discourse topic_id returned when the topic was created. Used when creating new posts for topics.

    """
    topic_title = models.CharField(max_length=200, unique=True,
                                   validators=
                                   [MinLengthValidator(15,'Discourse Topic Title must be longer than 15 characters')])
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    discourse_topic_id = models.IntegerField()

    class Meta:
        db_table = 'viewer_discoursetopic'
# End of Discourse Tables


class DownloadLinks(models.Model):
    """Django model containing the searches made with the download_structures
    api.

    Parameters
    ----------
    file_url: charField
        Contains the complete link to the zip file including the uuid.
    user: FK (integer)
        The (Django) id of the user that created the search
    target: FK (integer)
        The id of the target to which the tag belongs
    proteins: JSONField
        JSON field containing a sorted list of the protein codes in the search
    protein_params: JSONField
        JSON field containing sorted list of parameters used to create the
        zip file
    other_params: JSONField
        JSON field containing sorted list of parameters used to create the
        zip file
    static_link: BooleanField
        This preserves the proteins from the previous search.
    zip_contents: JSONField
        For static files, this field contains the contents of the zip so that
        it can be reconstructed with the same file-links that it had previously.
        For dynamic files, the zip is reconstructed from the search.
    create_date: DateTimeField
        The datetime when the search was created
    keep_zip_until: DateTimeField
        The datetime when the tag was created plus the retention time (1 hour
        at the time of writing)
    zip_file: BooleanField
        Link to the zip file created as part of the search - can be false if
        after the keep_zip_until.

    """
    file_url = models.CharField(max_length=200, unique=True, db_index=True)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE, db_index=True)
    proteins = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    protein_params = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    other_params = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    static_link = models.BooleanField(default=False)
    zip_contents = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    create_date = models.DateTimeField()
    keep_zip_until = models.DateTimeField(db_index=True)
    zip_file = models.BooleanField(default=False)
    original_search = models.JSONField(encoder=DjangoJSONEncoder, null=True)

    class Meta:
        db_table = 'viewer_downloadlinks'


# Start of Tag Tables
class TagCategory(models.Model):
    """Django model containing categories for Tags

    Parameters
    ----------
    category: CharField
        The unique name of the tag category.
    colour: CharField
        Expected to be an RGB string
    description: CharField
        Expected to be a helpful description of what will be contained in the tag category

    """
    category = models.CharField(max_length=50, unique=True)
    colour = models.CharField(max_length=20, null=True)
    description = models.CharField(max_length=200, null=True)

    class Meta:
        db_table = 'viewer_tagcategory'


class Tag(models.Model):
    """Django model containing Tags

    Parameters
    ----------
    tag: CharField
        The unique name of the tag.
    category: FK (integer)
        The id of the tag category to which the tag belongs
    target: FK (integer)
        The id of the target to which the tag belongs
    user: FK (integer)
        The (Django) id of the user that created the tag
    create_date: DateTimeField
        The datetime when the tag was created
    colour: CharField
        Expected to be an RGB string
    discourse_url: TextField
        Optional URL of a related Discourse Post
    help_text: TextField
        Optional help text to for the tag
    additional_info: JSONField
        Optional JSON field containing name/value pairs for future use

    """
    tag = models.CharField(max_length=200)
    category = models.ForeignKey(TagCategory, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    create_date = models.DateTimeField(default=timezone.now)
    colour = models.CharField(max_length=20, null=True)
    discourse_url = models.TextField(max_length=1000, null=True)
    help_text = models.TextField(null=True)
    additional_info = models.JSONField(encoder=DjangoJSONEncoder, null=True)

    class Meta:
        abstract = True
        unique_together = ('tag', 'target',)


class MoleculeTag(Tag):
    """Django model containing data for MoleculeTag(s) inherited from Tag.

    Parameters
    ----------
    molecules: ManyToManyField
        Links to the Molecule(s) that are tagged
    mol_group: ForeignKey scoring.Molgroup
        Links to the Molecule group - used for Sites when reloading Molecules

    """
    molecules = models.ManyToManyField(Molecule, blank=True)
    mol_group = models.ForeignKey("scoring.MolGroup", null=True, blank=True,
                                  on_delete=models.SET_NULL)


class SessionProjectTag(Tag):
    """Django model containing data for SessionProjectTag(s) inherited from Tag.

    Parameters
    ----------
    sesssion_peojects: ManyToManyField
        Links to the Session Projects) that are tagged

    """
    session_projects = models.ManyToManyField(SessionProject)
# End of Tag Tables


# Start of Squonk Job Tables
class JobFileTransfer(models.Model):
    """Django model containing Squonk File transfer contents and status
    information

    Parameters
    ----------
    id: Autofield
        Auto-created id for the file transfer.
    user: ForeignKey
        Foreign key link to the id of the user that created the file transfer request
    snapshot: ForeignKey
        A foreign key link to the relevant snapshot the file transfer is part of (required)
    target: ForeignKey
        A foreign key link to the relevant target the file transfer is part of (required)
    squonk_project: CharField
        The name of a project that has been created in Squonk that the files will be transferred to
    proteins: JSONField
        List of proteins to be transferred
    compounds: JSONField
        List of compounds to be transferred (not used yet)
    transfer_spec: JSONField
        Identifies for each type (protein or compound), which file types were transferred over.
    transfer_task_id: CharField
        Task id of transfer celery task Note that if a re-synchronisation is
        required this will be re-used.
    transfer_status: CharField
        Identifies the status of the transfer.
    transfer_progress: DecimalField
        Intended to be used as an indication of progress (0 to 100%)
    transfer_datetime: DateTimeField
        The datetime the transfer was completed.
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
    id = models.AutoField(primary_key=True)
    user = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE, db_index=True)
    squonk_project = models.CharField(max_length=200, null=True)
    sub_path = ShortUUIDField(length=4, alphabet="abcdefghijklmnopqrstuvwxyz")
    proteins = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    # Not used in phase 1
    compounds = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    transfer_spec = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    transfer_task_id = models.CharField(null=True, max_length=50)
    transfer_status = models.CharField(choices=STATUS, default=PENDING, max_length=7)
    transfer_progress = models.DecimalField(null=True, max_digits=5, decimal_places=2)
    transfer_datetime = models.DateTimeField(null=True)

    class Meta:
        db_table = 'viewer_jobfiletransfer'


class JobRequest(models.Model):
    """Django model containing Squonk Job Request information.
    Note that this will be updated from both the Fragalysis frontend
    (when the user will be logged on) and Squonk (using a call back URL)

    Parameters
    ----------
    id: Autofield
        Auto-created id for the file transfer.
    squonk_job_name: CharField
        The name of the squonk job that will be run
    user: ForeignKey
        Foreign key link to the id of the user that created the file transfer request
    snapshot: ForeignKey
        A foreign key link to the relevant snapshot the file transfer is part of (required)
    target: ForeignKey
        A foreign key link to the relevant target the file transfer is part of (required)
    project: ForeignKey
        The Fragalysis Project record ths JobRequest is associated with
    squonk_project: CharField
        The name of a project that has been created in Squonk that the files will be transferred to
    squonk_job_spec: JSONField
        The specification of the job that will be provided to the Squonk POST instance API
    job_start_datetime: DateField
        The datetime when the Squonk Job has started, populated by information in the
        Squonk callback.
    job_finish_datetime: DateField
        The datetime when the Squonk Job has finished, populated by information in the
        Squonk callback. If this is not set you can assume the JOb is still running.
        When it is set the job_status filed will be updated (to SUCCESS or FAILURE).
        If automatic upload follows an upload_task_id wil be set and you can monitor
        upload_status for a status of the upload
    job_status: CharField
        The status of the Squonk job. Will be modified by Squonk through the callback URL
    job_status_datetime: DateField
        The datetime of the most recent job_status change.
    squonk_job_info: JSONField
        Squonk job information returned from the initial Squonk POST instance API call
    squonk_url_ext: CharField
        Squonk URL information to be added to the Host URL to link to a Squonk Job.
        This field is populated during a call to the JobCallbackView.
    code: UUIDField
        A UUID generated by Fragalysis and passed to Squonk as part of a callback URL.
        This value is used to uniquely identify the HJob in Squonk and is passed back
        by squonk to provide context in calls to the JobCallBackView.
        context in subs
    upload_task_id: CharField
        Celery task ID for results upload task (optional). Set when the Job completes
        and an automated upload follows.
    upload_status: CharField
        Status for results upload task (optional)
    computed_set: ForeignKey
        ID of uploaded computed set (optional)
    """
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
    snapshot = models.ForeignKey(Snapshot, on_delete=models.CASCADE)
    target = models.ForeignKey(Target, null=True, on_delete=models.CASCADE,
                               db_index=True)
    project = models.ForeignKey(Project, null=True, on_delete=models.CASCADE)
    squonk_project = models.CharField(max_length=200, null=True)
    squonk_job_spec = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    # Start and finish times for the Job
    job_start_datetime = models.DateTimeField(null=True)
    job_finish_datetime = models.DateTimeField(null=True)
    # Job status (and status Datetime), delivered via callbacks from Squonk
    job_status = models.CharField(choices=SQUONK_STATUS, default=PENDING, max_length=7)
    job_status_datetime = models.DateTimeField(null=True)
    # squonk_job_info is a copy of the response from DmApi.start_job_instance().
    # It's an instance of a DmApiRv object (a namedtuple)
    # that contains a 'success' (boolean) and 'msg' (the DmApi response's resp.json()).
    # For us this will contain a 'task_id', 'instance_id' and 'callback_token'.
    # The content will be a list with index '0' that's the value of the DmApiRv
    # 'success' variable and, at index '1', the original response message json().
    # The Job callback token will be squonk_job_info[1]['callback_token']
    squonk_job_info = models.JSONField(encoder=DjangoJSONEncoder, null=True)
    # 'squonk_url_ext' is a Squonk UI URL to obtain information about the
    # running instance. It's essentially the Squonk URL with the instance ID appended.
    squonk_url_ext = models.CharField(max_length=200, null=True)
    code = models.UUIDField(default=uuid.uuid4, editable=False, unique=True)
    upload_task_id = models.CharField(null=True, max_length=50)
    upload_status = models.CharField(choices=UPLOAD_STATUS, default=PENDING, max_length=7,
                                     null=True)
    computed_set = models.ForeignKey(ComputedSet, on_delete=models.CASCADE, null=True)

    class Meta:
        db_table = 'viewer_jobrequest'

    def job_has_finished(self):
        """Finished if status is SUCCESS or FAILURE (or a new state of LOST).
        """
        return self.job_status in [JobRequest.SUCCESS, JobRequest.FAILURE, 'LOST']

class Squonk2Org(models.Model):
    """Django model to store Squonk2 Organisations (UUIDs) and the Account Servers
    they belong to. Managed by the Squonk2Agent class.

       Parameters
       ----------
       uuid: TextField (40)
           A Squonk2 Account Server (AS) Organisation UUID. A fixed length string
           consisting of 'org-' followed by a uuid4 value,
           e.g. 'org-54260047-183b-42e8-9658-385a1e1bd236'
       name: TextField (80)
           The name of the Squonk2 Organisation UUID (obtained form the AS).
       as_url: URLField (200)
           The URL of the Squonk2 Account Server that owns the organisation.
           e.g. 'https://example.com/account-server-api'
       as_version: TextField
           The version of the AS that was first seen to own the Organisation
       """
    uuid = models.TextField(max_length=40, null=False)
    name = models.TextField(max_length=80, null=False)
    as_url = models.URLField(null=False)
    as_version = models.TextField(null=False)

class Squonk2Unit(models.Model):
    """Django model to store Squonk2 Unit (UUIDs). Managed by the Squonk2Agent class.

       Parameters
       ----------
       uuid: TextField (41)
           A Squonk2 Account Server (AS) Unit UUID. A fixed length string
           consisting of 'unit-' followed by a uuid4 value,
           e.g. 'unit-54260047-183b-42e8-9658-385a1e1bd236'
       name: TextField (80)
           The name used to create the Squonk2 Unit UUID
           This is not limited by the actual name length imposed by the DM
       target_access: ForeignKey
           A Foreign Key to the Project (Proposal) the Unit belongs to,
           a record that contains the "target access string".
       organisation: ForeignKey
           The Organisation the Unit belongs to.
       """
    uuid = models.TextField(max_length=41, null=False)
    name = models.TextField(null=False)

#    target_access = models.ForeignKey(Project, null=False, on_delete=models.CASCADE)
    organisation = models.ForeignKey(Squonk2Org, null=False, on_delete=models.CASCADE)

class Squonk2Project(models.Model):
    """Django model to store Squonk2 Project and Product (UUIDs).
    Managed by the Squonk2Agent class.

       Parameters
       ----------
       uuid: TextField (44)
           A Squonk2 Data Manager (DM) Project UUID. A fixed length string
           consisting of 'project-' followed by a uuid4 value,
           e.g. 'project-54260047-183b-42e8-9658-385a1e1bd236'
       name: TextField (80)
           The name of the Squonk2 Unit UUID (obtained form the AS).
       product_uuid: TextField (44)
           A Squonk2 Account Server (AS) Product UUID. A fixed length string
           consisting of 'product-' followed by a uuid4 value,
           e.g. 'product-54260047-183b-42e8-9658-385a1e1bd236'
       unit: ForeignKey
           The Squonk2 Unit the Product (and Project) belongs to.
       """
    uuid = models.TextField(max_length=44, null=False)
    name = models.TextField(null=False)
    product_uuid = models.TextField(max_length=44, null=False)

    unit = models.ForeignKey(Squonk2Unit, null=False, on_delete=models.CASCADE)
#    user = models.ForeignKey(User, null=False, on_delete=models.CASCADE)
#    session_project = models.ForeignKey(SessionProject, null=False, on_delete=models.CASCADE)

# End of Squonk Job Tables
