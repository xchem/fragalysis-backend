from xchem_db.models import Crystal

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


class CrystalSerializer(serializers.ModelSerializer):

    class Meta:
        Model = Crystal
        fields = (
            'compound',
        )


