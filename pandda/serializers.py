from rest_framework import serializers
from pandda.models import PanddaSite, PanddaEvent

class PanddaEventSerializer(serializers.ModelSerializer):

    class Meta:
        model = PanddaEvent
        fields = ('id','xtal','event','pandda_site','target_id','pdb_info','mtz_info','map_info','lig_id',
                  'event_com_x','event_com_y','event_com_z','lig_com_x','lig_com_y','lig_com_z',
                  'event_dist_from_site_centroid','lig_dist_from_site_centroid','small_map_info')

class PanddaSiteSerializer(serializers.ModelSerializer):

    class Meta:
        model = PanddaSite
        fields = ('id','site_id','pandda_run','pandda_version','target_id','site_align_com_x','site_align_com_y','site_align_com_z',
                  'site_native_com_x','site_native_com_y','site_native_com_z',)