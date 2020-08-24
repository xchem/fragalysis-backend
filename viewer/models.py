from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone
import uuid
import json

from simple_history.models import HistoricalRecords

from loader.config import get_mol_choices, get_prot_choices


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
    """
    # The title of the project_id -> userdefined
    title = models.CharField(max_length=200, unique=True)
    # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    # The users it's related to
    user_id = models.ManyToManyField(User)


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
    """
    # The title of the project_id -> userdefined
    title = models.CharField(unique=True, max_length=200)
    # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    # A field to link projects and targets together
    project_id = models.ManyToManyField(Project)
    # Indicates the uniprot_id id for the target. Is a unique key
    uniprot_id = models.CharField(max_length=100, null=True)
    # metadatafile containing sites info for download
    metadata = models.FileField(upload_to="metadata/", null=True, max_length=255)
    # zip archive to download uploaded data from
    zip_archive = models.FileField(upload_to="archive/", null=True, max_length=255)


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
    aligned: NullBooleanField
        Bool - 1 if aligned, 0 if not
    aligned_to: ForeignKey (self)
        Foreign key to another instance of a protein to which this protein is aligned (optional)
    has_eds: NullBooleanField
        Bool - 1 if has ED, 0 it not (optional)
    """
    # code for this protein (e.g. NUDT5A-x0001_1)
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
    aligned = models.NullBooleanField()
    aligned_to = models.ForeignKey("self", null=True, on_delete=models.CASCADE)
    has_eds = models.NullBooleanField()

    class Meta:
        unique_together = ("code", "prot_type")


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
    inchi = models.CharField(max_length=255, unique=True, db_index=True)
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


# Start of Session Project
class SessionProject(models.Model):
    """Django model for holding information about a fragalysis user project - a set of sessions saved by a user

    Parameters
    ----------
    title: CharField
        The title of the project
    init_date: DateTimeField
        The date the project was initiated (autofield)
    description: Charfield
        A short user-defined description for the project
    target: ForeignKey
        Foreign Key link to the relevent project target
    author: ForeignKey
        A link to the user that created the project
    tags: TextField
        A comma separated list of user-defined tags - for searching and tagging projects

    """
    title = models.CharField(max_length=200)
    init_date = models.DateTimeField(default=timezone.now)
    description = models.CharField(max_length=255, default='')
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    author = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    tags = models.TextField(default='[]')

    class Meta:
        db_table = 'viewer_sessionproject'


class Snapshot(models.Model):
    """Django model for storing information describing a static snapshot of a fragalysis page - multiple snapshots make
    up a project

    Parameters
    ----------
    id: Autofied
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

    class Meta:
        managed = True
        db_table = 'viewer_snapshot'


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
    """
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
    # the protein link
    pdb_info = models.FileField(upload_to="pdbs/", null=False, max_length=255)
    # if we use our own method of calculating them
    computed_inspirations = models.ManyToManyField(Molecule, null=True, blank=True)


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
