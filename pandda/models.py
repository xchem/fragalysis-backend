from django.db import models
from viewer.models import Target

class PanddaSite(models.Model):
    """
    Django model for Pandda sites
    """
    # The site
    site_id = models.IntegerField()
    # The pdb, map and mtz
    target_id = models.ForeignKey(Target)
    # The pannda version
    pandda_version = models.TextField(null=True)
    # Pandda run
    pandda_run = models.TextField(default="STANDARD")
    # Define the X,Y,Z
    site_align_com_x = models.FloatField(null=True)
    site_align_com_y = models.FloatField(null=True)
    site_align_com_z = models.FloatField(null=True)
    site_native_com_x = models.FloatField(null=True)
    site_native_com_y = models.FloatField(null=True)
    site_native_com_z = models.FloatField(null=True)

    class Meta:
        unique_together = ('site_id', 'target_id', 'pandda_run')

class PanddaEvent(models.Model):
    # Define the site
    pandda_site = models.ForeignKey(PanddaSite)
    target_id = models.ForeignKey(Target)
    # The xtal - this will later be linke to Rachael's stuff
    xtal = models.TextField()
    # The event id
    event = models.IntegerField()
    apo_holo = models.NullBooleanField()
    pdb_info = models.FileField(upload_to='pdbs/', null=True, max_length=10000000)
    cif_info = models.FileField(upload_to='cifs/', null=True, max_length=10000000)
    mtz_info = models.FileField(upload_to='mtzs/', null=True, max_length=10000000)
    map_info = models.FileField(upload_to='maps/', null=True, max_length=10000000)
    small_map_info = models.FileField(upload_to='maps/', null=True, max_length=10000000)
    # The ligand id
    lig_id = models.CharField(max_length=50, null=True)
    event_com_x = models.FloatField(null=True)
    event_com_y = models.FloatField(null=True)
    event_com_z = models.FloatField(null=True)
    lig_com_x = models.FloatField(null=True)
    lig_com_y = models.FloatField(null=True)
    lig_com_z = models.FloatField(null=True)
    event_dist_from_site_centroid = models.FloatField(null=True)
    lig_dist_from_site_centroid = models.FloatField(null=True)
    # Unique constraints
    class Meta:
        unique_together = ('xtal', 'event', 'pandda_site', 'target_id')