from django.db import models
from django.contrib.auth.models import User
from viewer.models import Target,Protein,Molecule
from analysis.models import Fragment, Pharmacaphore

class ResShift(models.Model):
    """Model to store residue shift information"""
    # The target it relates to
    target_id = models.ForeignKey(Target)
    # The residue number
    res_num = models.IntegerField()
    # The residue name
    res_name = models.CharField(max_length=20)
    # The max shift for this residue
    max_shift = models.FloatField(null=True)
    # The average shift for this residue
    avg_shift = models.FloatField(null=True)
    # The min shift for this residue
    min_shift = models.FloatField(null=True)

    class Meta:
        unique_together = ('target_id', 'res_num', 'res_name')


class Residue(models.Model):
    """Model to store residue information - to curate the probes"""
    # The target it relates to
    target_id = models.ForeignKey(Target)
    # The protein it relates to
    prot_id = models.ForeignKey(Protein)
    # The residue number
    res_num = models.IntegerField()
    # The residue name
    res_name = models.CharField(max_length=20)
    # The x_min, x_max, y_min, y_max and z_min, z_max
    x_min = models.FloatField(null=True)
    x_max = models.FloatField(null=True)
    y_min = models.FloatField(null=True)
    y_max = models.FloatField(null=True)
    z_min = models.FloatField(null=True)
    z_max = models.FloatField(null=True)
    # The max shift for this residue
    max_shift = models.FloatField(null=True)
    # The average shift for this residue
    avg_shift = models.FloatField(null=True)
    # The min shift for this residue
    min_shift = models.FloatField(null=True)
    clust_id = models.IntegerField(null=True)
    # The PDB informaton  of this residue
    pdb_info = models.TextField()

    class Meta:
        unique_together = ('prot_id', 'res_num', 'res_name')


class InternalIDLink(models.Model):
    """Model to link a molecule to an internal ID and therefore to activity data"""
    # The molecule it links to
    mol_id = models.ForeignKey(Molecule)
    # The ID of the compound for internal use - to reference the activity data
    internal_id = models.CharField(max_length=150, null=True)


class PharmaScore(models.Model):
    """A score for a pharmapoint"""
    score = models.FloatField()
    type = models.CharField(max_length=30)
    pharma_id = models.ForeignKey(Pharmacaphore)

    class Meta:
        unique_together = ('pharma_id', 'type', )


class FragScore(models.Model):
    """A score for an MMPFrag"""
    score = models.FloatField()
    type = models.CharField(max_length=30)
    frag_id = models.ForeignKey(Fragment)

    class Meta:
        unique_together = ('frag_id', 'type', )


class Water(models.Model):
    """A water molecule"""
    # The coords
    x_com = models.FloatField(db_index=True)
    y_com = models.FloatField(db_index=True)
    z_com = models.FloatField(db_index=True)
    # The target it corresponds to
    target_id = models.ForeignKey(Target)
    # The protein it refers to
    prot_id = models.ForeignKey(Protein)
    # The identifier for that water
    water_num = models.IntegerField()

    class Meta:
        unique_together = ('prot_id', 'water_num', 'target_id')


class KeyCluster(models.Model):
    """A key cluster for considering"""
    # The target it corresponds to
    target_id = models.ForeignKey(Target)
    # The number of features associated
    size = models.IntegerField(null=True, db_index=True)
    # Is it a ph4 or a frag
    type = models.CharField(max_length=300, db_index=True)
    # Is it an acceptor or whatever
    function = models.CharField(max_length=300, db_index=True)
    # The lamda used
    lam = models.FloatField(db_index=True)
    # The id of this cluster
    cluster = models.IntegerField(db_index=True)
    # The method used to cluster
    method = models.CharField(max_length=300, db_index=True)
    # Boolean to show a fragment is just a ring system
    isring = models.BooleanField(default=False)
    frag_size = models.IntegerField(default=0, db_index=True)
    # The coords
    x_com = models.FloatField(db_index=True)
    y_com = models.FloatField(db_index=True)
    z_com = models.FloatField(db_index=True)
    # The pharmapoint associated to it
    pharma_id = models.ManyToManyField(Pharmacaphore)
    # The MMPFrag associate to it
    mmp_frag_id = models.ManyToManyField(Fragment)
    # The water associated to it
    water_id = models.ManyToManyField(Water)
    # The molecules associated to it
    mol_id = models.ManyToManyField(Molecule)
    # The cluster max and min
    x_max = models.FloatField(db_index=True, null=True)
    x_min = models.FloatField(db_index=True, null=True)
    y_max = models.FloatField(db_index=True, null=True)
    y_min = models.FloatField(db_index=True, null=True)
    z_max = models.FloatField(db_index=True, null=True)
    z_min = models.FloatField(db_index=True, null=True)
    # Ringname
    ringname = models.TextField(db_index=True, null=True)

    class Meta:
        unique_together = ('type', 'lam', 'method', 'function', 'cluster', 'target_id')
