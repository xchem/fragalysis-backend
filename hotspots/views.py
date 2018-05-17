from hotspots.models import HotspotMap
from hotspots.serializers import HotspotMapSerializer
from rest_framework import permissions
from rest_framework import viewsets

class HotspotView(viewsets.ReadOnlyModelViewSet):
    queryset = HotspotMap.objects.filter()
    serializer_class = HotspotMapSerializer
    filter_fields =  ('map_type','target_id','prot_id',)