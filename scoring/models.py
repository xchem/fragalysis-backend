from django.contrib.auth.models import User
from django.db import models

from viewer.models import SiteObservation, Compound, Target, Snapshot


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


class SiteObservationChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
    DEFAULT = "DE"
    PANDDA = "PA"
    GOOD_MOL = "GM"
    # TODO: what choices will there be?
    MOL_CHOICES = (
        (DEFAULT, "Default"),
        (PANDDA, "Pandda"),
        (GOOD_MOL, "Good molecule"),
    )
    choice_type = models.CharField(choices=MOL_CHOICES, max_length=2, default=DEFAULT)
    # Score -
    score = models.FloatField(null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["user_id", "site_observation", "choice_type"],
                name="unique_siteobvs_choice",
            ),
        ]         
        

class SiteObservationAnnotation(models.Model):
    """
    A Django model to annotate a molecule with free text
    """
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
    annotation_type = models.CharField(max_length=50)
    annotation_text = models.CharField(max_length=100)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["site_observation", "annotation_type"],
                name="unique_siteobvs_anntype",
            ),
        ]         


class ScoreChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    # IN THIS CASE THIS WOULD INDICATE THE SOFTWARE USED - WE WILL GENERATE DIFFERENT USERS FOR EACH SOFTWARE
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
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
        constraints = [
            models.UniqueConstraint(
                fields=["user", "site_observation", "choice_type"],
                name="unique_user_siteobvs_dockchoice",
            ),
        ]        


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


class SiteObservationGroup(models.Model):
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
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
    # Set the description
    description = models.TextField(null=True)
    # Set the ManyToMany
    site_observation = models.ManyToManyField(
        SiteObservation,
        through="SiteObvsSiteObservationGroup",
        through_fields=("site_obvs_group", "site_observation"),
        related_name="site_observation_groups",
    )
    # Set the centre of mass
    x_com = models.FloatField(null=True)
    y_com = models.FloatField(null=True)
    z_com = models.FloatField(null=True)

    def __str__(self) -> str:
        return f"{self.pk}"

    def __repr__(self) -> str:
        return "<Xtalform %r %r %r>" % (self.id, self.group_type, self.description)


class SiteObvsSiteObservationGroup(models.Model):
    site_obvs_group = models.ForeignKey(SiteObservationGroup, null=False, on_delete=models.CASCADE)
    site_observation = models.ForeignKey(SiteObservation, null=False, on_delete=models.CASCADE)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["site_observation", "site_obvs_group",],
                name="unique_siteobservationgroupcontents",
            ),
        ]    
    
