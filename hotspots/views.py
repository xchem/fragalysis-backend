from rest_framework import viewsets

from hotspots.models import HotspotMap
from hotspots.serializers import HotspotMapSerializer


class HotspotView(viewsets.ReadOnlyModelViewSet):
    paginate_by = None
    queryset = HotspotMap.objects.filter()
    serializer_class = HotspotMapSerializer
    filter_fields = ("map_type", "target_id", "prot_id")
