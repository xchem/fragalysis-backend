from django.db import models
from django.contrib.auth.models import User
from viewer.models import Protein,Molecule,Compound


# TODO - Consider this for reorginisation. E.g. Linking to targets,proteins,molecules.
class ViewScene(models.Model):
    """
    A Django Model for a given view.
    This probably needs cleaning up - linking

    """
    # The user who made the scene
    user_id = models.ForeignKey(User,null=True)
    # The uuid that will enable this to be shared
    uuid = models.UUIDField(unique=True)
    # The title - optional
    title = models.CharField(max_length=200,default="NA")
    # The JSON describing the data
    scene = models.TextField()

    permissions = (
        ('view_scene', 'View scene'),
    )

class ProtChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user_id = models.ForeignKey(User)
    prot_id = models.ForeignKey(Protein)
    # Set the groups types
    DEFAULT = 'DE'
    PROT_CHOICES = (
        (DEFAULT, 'Default'),
    )
    choice_type = models.CharField(choices=PROT_CHOICES,max_length=2,default=DEFAULT)
    # Integer Score for this
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ('user_id', 'prot_id', 'choice_type',)
        permissions = (
            ('view_protchoice', 'View protchoice'),
        )


class MolChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user_id = models.ForeignKey(User)
    mol_id = models.ForeignKey(Molecule)
    DEFAULT = 'DE'
    PANDDA = 'PA'
    GOOD_MOL = 'GM'
    MOL_CHOICES = (
        (DEFAULT, 'Default'),
        (PANDDA, 'Pandda'),
        (GOOD_MOL, 'Good molecule'),
    )
    choice_type = models.CharField(choices=MOL_CHOICES,max_length=2,default=DEFAULT)
    # Score -
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ('user_id', 'mol_id', 'choice_type')
        permissions = (
            ('view_molchoice', 'View molchoice'),
        )


class ScoreChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    # IN THIS CASE THIS WOULD INDICATE THE SOFTWARE USED - WE WILL GENERATE DIFFERENT USERS FOR EACH SOFTWARE
    user_id = models.ForeignKey(User)
    mol_id = models.ForeignKey(Molecule)
    prot_id = models.ForeignKey(Molecule)
    is_done = models.BooleanField(default=False)
    DEFAULT = 'DE'
    DOCKING = 'AU'
    DENSITY_FIT = 'DF'
    INTERACTION = "IT"
    DOCK_CHOICES = (
        (DEFAULT, 'Default'),
        (DOCKING, 'Docking'),
        (INTERACTION, 'Interaction fit'),
        (DENSITY_FIT, 'Density Fit'),
    )
    choice_type = models.CharField(choices=DOCK_CHOICES,max_length=2,default=DEFAULT)
    # Any score
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ('user_id', 'mol_id', 'prot_id', 'choice_type')
        permissions = (
            ('view_scorechoice', 'View scorechoice'),
        )

class CmpdChoice(models.Model):
    """
    A Django model to store a selection from a user
    """
    user_id = models.ForeignKey(User)
    cmpd_id = models.ForeignKey(Compound)
    DEFAULT = 'DE'
    PRICE = 'PR'
    TOXIC = 'TO'
    CMPD_CHOICES = (
        (DEFAULT, 'Default'),
        (PRICE, 'Price'),
        (TOXIC, 'Toxic'),
    )
    choice_type = models.CharField(choices=CMPD_CHOICES,max_length=2,default=DEFAULT)
    # Score between 0 and 9; Convention for n memberd list -> num_in_list/(num_choices-1)
    # E.g.
    score = models.FloatField(null=True)

    class Meta:
        unique_together = ('user_id', 'cmpd_id','choice_type',)
        permissions = (
            ('view_cmpdchoice', 'View cmpdhoice'),
        )