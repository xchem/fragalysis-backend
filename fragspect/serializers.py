from xchem_db.models import Crystal, PanddaSite, PanddaEvent, PanddaEventStats, DataProcessing

from rest_framework import serializers

# "fragId": 55,
# "crystal": "NUDT7A_Crude-x0526",
# "site_number": "1",
# "event_number": "1",
# "code": "NUDT7A_Crude-x0526_3",
# "lig_id": "LIG-E1",
# "target_name": "NUDT7A",
# "target_id": 5,
# "prot_id": 8655,
# "event_map_info": "media/maps/MURD-x0349_event_ebBZqDc.ccp4",
# "sigmaa_map_info": "media/maps/MURD-x0349_sigmaa_ebBZqDc.ccp4",
# "spider_plot_info": "media/spideys/MURD-x0349_spider_ebBZqDc.png",
# "2d_view_png": "media/2dv/MURD-x0349_2Dview_ebBZqDc.png",
# "pandda_model_found": true,
# "crystal_status": 6,
# "event_status": 2,
# "confidence": 3,
# "resolution": 2.5,
# "smiles": "Cc1cc(NC(=O)Cc2cccc(O)c2)no1",
# "space_group": "C 1 2 1",
# "cell_dimensions": "102, 45, 60",
# "cell_angles": "90, 90, 90"


class FragspectSerializer(serializers.ModelSerializer):
    site_number = serializers.SerializerMethodField()
    event_number = serializers.SerializerMethodField()
    # codes = serializers.SerializerMethodField()
    lig_id = serializers.SerializerMethodField()
    target_name = serializers.SerializerMethodField()
    smiles = serializers.SerializerMethodField()
    # pandda_model_found = serializers.SerializerMethodField()
    # resolution = serializers.SerializerMethodField()
    space_group = serializers.SerializerMethodField()
    unit_cell = serializers.SerializerMethodField()

    def get_site_number(self, obj):
        return [e.site.site for e in PanddaEvent.objects.filter(crystal=obj)]

    def get_event_number(self, obj):
        return [e.event for e in PanddaEvent.objects.filter(crystal=obj)]

    # def get_code(self, obj):
    #     pass

    def get_lig_id(self, obj):
        return [e.lig_id for e in PanddaEvent.objects.filter(crystal=obj)]

    def get_target_name(self, obj):
        return [e.crystal.target.target_name for e in PanddaEvent.objects.filter(crystal=obj)]

    def get_smiles(self, obj):
        return [e.crystal.compound.smiles for e in PanddaEvent.objects.filter(crystal=obj)]

    # def get_pandda_model_found(self, obj):
    #     pass

    # def get_resolution(self, obj):
    #     try:
    #         return [PanddaEventStats.objects.get(event=e).resolution for e in PanddaEvent.objects.filter(crystal=obj)]
    #     except:


    def get_space_group(self, obj):
        return [DataProcessing.objects.get(crystal_name=e.crystal).spacegroup
                for e in PanddaEvent.objects.filter(crystal=obj)]

    def get_unit_cell(self, obj):
        return [DataProcessing.objects.get(crystal_name=e.crystal).unit_cell
                for e in PanddaEvent.objects.filter(crystal=obj)]

    class Meta:
        Model = Crystal
        fields = (
            'site_number',
            'event_number',
            'lig_id',
            'target_name',
            'smiles',
            'resolution',
            'spacegroup',
            'unit_cell',
        )


