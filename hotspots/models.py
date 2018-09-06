from django.db import models

from viewer.models import Target, Protein


class HotspotMap(models.Model):
    """
    Django model for Hotspot Maps
    """
    # The site
    prot_id = models.ForeignKey(Protein)
    # The pdb, map and mtz
    target_id = models.ForeignKey(Target)
    # The type of map
    APOLAR = "AP"
    ACCEPTOR = "AC"
    DONOR = "DO"
    MAP_CHOICES = ((ACCEPTOR, "Acceptor"), (APOLAR, "Apolar"), (DONOR, "Donor"))
    # Set the groups types
    map_type = models.CharField(choices=MAP_CHOICES, max_length=2)
    # The map file itself
    map_info = models.FileField(upload_to="maps/", null=True, max_length=255)
    compressed_map_info = models.FileField(upload_to="maps/", null=True, max_length=255)

    class Meta:
        unique_together = ("map_type", "target_id", "prot_id")
