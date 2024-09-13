from api.security import ISPyBSafeQuerySet
from hotspots.models import HotspotMap
from hotspots.serializers import HotspotMapSerializer


class HotspotView(ISPyBSafeQuerySet):
    queryset = HotspotMap.objects.all()
    serializer_class = HotspotMapSerializer
    filterset_fields = ("map_type", "target", "site_observation")
    filter_permissions = "target__project_id"
