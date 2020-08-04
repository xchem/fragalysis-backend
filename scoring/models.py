from django.contrib.auth.models import User
from django.db import models

from viewer.models import Protein, Molecule, Compound, Target, Snapshot


class ViewScene(models.Model):
    """
    A Django Model for a given view.
    This probably needs cleaning up - linking

    """
    # The user who made the scene
    user_id = models.ForeignKey(User, null=True, on_delete=models.CASCADE)
    # The uuid that will enable this to be shared
    uuid = models.UUIDField(unique=True)
    # The title - optional
    title = models.CharField(max_length=200, default="NA")
    # The JSON describing the data
    scene = models.TextField()
    # autofield for when the scene was created
    created = models.DateTimeField(auto_now_add=True)
    # autofield for when the scene was modified
    modified = models.DateTimeField(auto_now=True)

    # for redirecting to project's snapshot
    snapshot = models.ForeignKey(Snapshot, null=True, on_delete=models.CASCADE)


class ProtChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user_id = models.ForeignKey(User, on_delete=models.CASCADE)
    prot_id = models.ForeignKey(Protein, on_delete=models.CASCADE)
    # Set the groups types
    DEFAULT = "DE"
    PROT_CHOICES = ((DEFAULT, "Default"),)
    choice_type = models.CharField(choices=PROT_CHOICES, max_length=2, default=DEFAULT)
    # Integer Score for this
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ("user_id", "prot_id", "choice_type")


class MolChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user_id = models.ForeignKey(User, on_delete=models.CASCADE)
    mol_id = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    DEFAULT = "DE"
    PANDDA = "PA"
    GOOD_MOL = "GM"
    MOL_CHOICES = (
        (DEFAULT, "Default"),
        (PANDDA, "Pandda"),
        (GOOD_MOL, "Good molecule"),
    )
    choice_type = models.CharField(choices=MOL_CHOICES, max_length=2, default=DEFAULT)
    # Score -
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ("user_id", "mol_id", "choice_type")


class MolAnnotation(models.Model):
    """
    A Django model to annotate a molecule with free text
    """
    mol_id = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    annotation_type = models.CharField(max_length=50)
    annotation_text = models.CharField(max_length=100)

    class Meta:
        unique_together = ("mol_id", "annotation_type")


class ScoreChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    # IN THIS CASE THIS WOULD INDICATE THE SOFTWARE USED - WE WILL GENERATE DIFFERENT USERS FOR EACH SOFTWARE
    user_id = models.ForeignKey(User, on_delete=models.CASCADE)
    mol_id = models.ForeignKey(Molecule, on_delete=models.CASCADE)
    prot_id = models.ForeignKey(Protein, on_delete=models.CASCADE)
    is_done = models.BooleanField(default=False)
    DEFAULT = "DE"
    DOCKING = "AU"
    DENSITY_FIT = "DF"
    INTERACTION = "IT"
    DOCK_CHOICES = (
        (DEFAULT, "Default"),
        (DOCKING, "Docking"),
        (INTERACTION, "Interaction fit"),
        (DENSITY_FIT, "Density Fit"),
    )
    choice_type = models.CharField(choices=DOCK_CHOICES, max_length=2, default=DEFAULT)
    # Any score
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ("user_id", "mol_id", "prot_id", "choice_type")


class CmpdChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user_id = models.ForeignKey(User, on_delete=models.CASCADE)
    cmpd_id = models.ForeignKey(Compound, on_delete=models.CASCADE)
    DEFAULT = "DE"
    PRICE = "PR"
    TOXIC = "TO"
    CMPD_CHOICES = ((DEFAULT, "Default"), (PRICE, "Price"), (TOXIC, "Toxic"))
    choice_type = models.CharField(choices=CMPD_CHOICES, max_length=2, default=DEFAULT)
    # Score between 0 and 9; Convention for n memberd list -> num_in_list/(num_choices-1)
    # E.g.
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ("user_id", "cmpd_id", "choice_type")


class MolGroup(models.Model):
    """
    A Django model for a group of molecules.
    No unique set - so needs to be deleted for a type before re-running.
    """
    PANDDA = "PA"
    DEFAULT = "DE"
    MOL_CLUSTER = "MC"
    WATER_CLUSTER = "WC"
    PHARMA_CLUSTER = "PC"
    RES_CLUSTER = "RC"
    MOL_GROUP_CHOICES = (
        (PANDDA, "Pandda"),
        (DEFAULT, "Default"),
        (MOL_CLUSTER, "MolCluster"),
        (WATER_CLUSTER, "WaterCluster"),
        (PHARMA_CLUSTER, "PharmaCluster"),
        (RES_CLUSTER, "ResCluster"),
    )
    # Set the groups types
    group_type = models.CharField(
        choices=MOL_GROUP_CHOICES, max_length=2, default=DEFAULT
    )
    # Set the target id
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    # Set the description
    description = models.TextField(null=True)
    # Set the ManyToMany
    mol_id = models.ManyToManyField(Molecule, related_name="mol_groups")
    # Set the centre of mass
    x_com = models.FloatField(null=True)
    y_com = models.FloatField(null=True)
    z_com = models.FloatField(null=True)

