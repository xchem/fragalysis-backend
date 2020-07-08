from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone
import uuid
import json

from simple_history.models import HistoricalRecords

from loader.config import get_mol_choices, get_prot_choices


class Project(models.Model):
    """A django model to define a given project_id. Not currently used.
    Could be used to define certain attributes."""
    # The title of the project_id -> userdefined
    title = models.CharField(max_length=200, unique=True)
    # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    # The users it's related to
    user_id = models.ManyToManyField(User)

    class Meta:
        permissions = (("view_project", "View project"),)


class Target(models.Model):
    """A Django model to define a given protein target"""
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

    class Meta:
        permissions = (("view_target", "View target"),)


class Protein(models.Model):
    """A Django model to hold the information for a given protein, unique set of coords"""
    # code for this protein (e.g. NUDT5A-x0001_1)
    code = models.CharField(max_length=50, db_index=True)
    target_id = models.ForeignKey(Target)
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
    aligned_to = models.ForeignKey("self", null=True)
    has_eds = models.NullBooleanField()

    class Meta:
        unique_together = ("code", "prot_type")
        permissions = (("view_protein", "View protein"),)


class Compound(models.Model):
    """A Django model to hold the information for a given compound -> a unique 2D molecule"""
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

    def set_all_identifiers(self, x):
        self.all_identifiers = json.dumps(x)

    def get_all_identifiers(self):
        return json.loads(self.foo)

    class Meta:
        permissions = (("view_compound", "View compound"),)
        # unique_together = (("inchi", "long_inchi"),)


class Molecule(models.Model):
    """A Django model to hold the information for a molecule -> a 3D set of
    co-ordinates"""
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
    prot_id = models.ForeignKey(Protein)
    cmpd_id = models.ForeignKey(Compound)
    history = HistoricalRecords()

    # Unique constraints
    class Meta:
        unique_together = ("prot_id", "cmpd_id", "mol_type")
        permissions = (("view_molecule", "View molecule"),)


class ActivityPoint(models.Model):
    """A Django model to hold the activity information for a given compound"""
    # This should encompass the activity type too
    source = models.CharField(max_length=50, null=True, db_index=True)
    # Foreign key to Target object
    target_id = models.ForeignKey(Target)
    # Foreign key to compound object
    cmpd_id = models.ForeignKey(Compound)
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
        permissions = (("view_activitypoint", "View activitypoint"),)



# Start of Session Project
class SessionProject(models.Model):
    title = models.CharField(max_length=200)
    init_date = models.DateTimeField(default=timezone.now)
    description = models.CharField(max_length=255, default='')
    target = models.ForeignKey(Target)
    author = models.ForeignKey(User, null=True)
    tags = models.TextField(default='[]')

    class Meta:
        permissions = (("view_project", "View project"),)
        db_table = 'viewer_sessionproject'

class Snapshot(models.Model):
    INIT = "INIT"
    AUTO = "AUTO"
    MANUAL = "MANUAL"
    SNAPSHOT_TYPE = (
        (INIT, "INIT"), # Initial snapshot generated by system
        (AUTO, 'AUTO'),  # Automatic generated by system
        (MANUAL,'MANUAL') # Manual generated by user action
    )
    id = models.AutoField(primary_key=True)
    type = models.CharField(choices=SNAPSHOT_TYPE, default=INIT, max_length=8)
    title = models.CharField(max_length=255)
    author = models.ForeignKey(User, null=True)
    description = models.CharField(max_length=255, default='')
    created = models.DateTimeField(default=timezone.now)
    data = models.TextField()
    session_project = models.ForeignKey(SessionProject, null=True)
    parent = models.ForeignKey('self', models.DO_NOTHING, blank=True, null=True, related_name='children')

    class Meta:
        permissions = (("view_project", "View project"),)
        managed = True
        db_table = 'viewer_snapshot'
# End of Session Project


# Start of compound sets
# Design sets = 2D compounds that have been designed but do not yet have a 3D structure
class DesignSet(models.Model):
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
    name = models.CharField(max_length=50, null=False)
    email = models.CharField(max_length=100, null=False)
    institution = models.CharField(max_length=50, null=False)
    generation_date = models.DateField()
    method = models.CharField(max_length=50, null=False)

    class Meta:
        unique_together = (("name", "method"),)


class CSetKeys(models.Model):
    user = models.CharField(max_length=50, default='User', editable=False)
    uuid = models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True)


# computed sets = sets of poses calculated computationally
class ComputedSet(models.Model):
    # a (unique) name for this compound set
    name = models.CharField(max_length=50, unique=True, primary_key=True)
    # target that this compound set belongs to
    target = models.ForeignKey(Target, null=True)
    # link to the submitted sdf file
    submitted_sdf = models.FileField(upload_to="compound_sets/", null=False, max_length=255)
    # file format specification version
    spec_version = models.FloatField(null=False)
    method_url = models.TextField(max_length=1000, null=True)
    submitter = models.ForeignKey(ComputedSetSubmitter, null=True)
    unique_name = models.CharField(max_length=101, null=False)
    # Check if needed? Rachael still having a look.
    #design_set = models.ForeignKey(DesignSet, null=False, blank=False)

    def save(self):
        if not self.submitter:
            super(ComputedSet, self).save()
        if not self.unique_name:
            unique_name = "".join(self.submitter.name.split()) + '-' + "".join(self.submitter.method.split())
            self.unique_name = unique_name
        super(ComputedSet, self).save()


class ComputedMolecule(models.Model):
    compound = models.ForeignKey(Compound)
    # the 3D coordinates of a computed molecule
    sdf_info = models.TextField(null=False)
    # a link to the compound set this molecule belongs to
    computed_set = models.ForeignKey(ComputedSet, on_delete=models.CASCADE)
    # a name for this compound
    name = models.CharField(max_length=50)
    # calculated smiles
    smiles = models.CharField(max_length=255)
    # comes from compound
    #original_smiles = models.CharField(max_length=255)
    # link to pdb file for prot structure
    pdb_info = models.FileField(upload_to="pdbs/", null=False, max_length=255)
    #design_set = models.ForeignKey(DesignSet) # needs to be linked to find inspiration fragments
    computed_inspirations = models.ManyToManyField(Molecule, null=True, blank=True) # if we use our own method of calculating them


class ScoreDescription(models.Model):
    # which compound set this score belongs to
    computed_set = models.ForeignKey(ComputedSet, null=True)
    # a name for this score
    name = models.CharField(max_length=50)
    # a description for this score
    description = models.TextField(null=False)


class NumericalScoreValues(models.Model):
    # a link to the score name and description
    score = models.ForeignKey(ScoreDescription)
    # the specific value for this score
    value = models.FloatField(null=False)
    # the compound this score relates to
    compound = models.ForeignKey(ComputedMolecule)


class TextScoreValues(models.Model):
    # a link to the score name and description
    score = models.ForeignKey(ScoreDescription)
    # the specific value for this score
    value = models.TextField(max_length=500, null=False)
    # the compound this score relates to
    compound = models.ForeignKey(ComputedMolecule)

# End of compound sets
class File(models.Model):
    file = models.FileField(blank=False, null=False)
    def __str__(self):
        return self.file.name


