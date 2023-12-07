from django.db import models

from viewer.models import SiteObservation, Target


class HotspotMap(models.Model):
    """
    Django model for Hotspot Maps
    """

    # The site
    # prot_id = models.ForeignKey(Protein, on_delete=models.CASCADE)
    site_observation = models.ForeignKey(SiteObservation, on_delete=models.CASCADE)
    # The pdb, map and mtz
    target = models.ForeignKey(Target, on_delete=models.CASCADE)
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
        constraints = [
            models.UniqueConstraint(
                fields=[
                    "map_type",
                    "target",
                    "site_observation",
                ],
                name="unique_maptype_target_siteobvs",
            ),
        ]
