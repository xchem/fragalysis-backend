from rest_framework import serializers

from hotspots.models import HotspotMap


class HotspotMapSerializer(serializers.ModelSerializer):
    class Meta:
        model = HotspotMap
        fields = (
            "id",
            "map_type",
            "target",
            "site_observation",
            "map_info",
            "compressed_map_info",
        )
