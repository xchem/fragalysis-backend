from rest_framework import viewsets

from hotspots.models import HotspotMap
from hotspots.serializers import HotspotMapSerializer


class HotspotView(viewsets.ReadOnlyModelViewSet):
    queryset = HotspotMap.objects.all()
    serializer_class = HotspotMapSerializer
    filterset_fields = ("map_type", "target", "site_observation")
