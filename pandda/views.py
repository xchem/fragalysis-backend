from pandda.models import PanddaEvent,PanddaSite
from pandda.serializers import PanddaEventSerializer,PanddaSiteSerializer
from rest_framework import permissions
from rest_framework import viewsets

class PanddaEventView(viewsets.ReadOnlyModelViewSet):
    queryset = PanddaEvent.objects.filter()
    serializer_class = PanddaEventSerializer
    filter_fields = ('xtal','event','pandda_site','target_id',)

class PanddaSiteView(viewsets.ReadOnlyModelViewSet):
    queryset = PanddaSite.objects.filter()
    serializer_class = PanddaSiteSerializer
    filter_fields = ('pandda_run','site_id','pandda_version','target_id',)