from xchem_db.models import Crystal

from rest_framework import serializers


class CrystalSerializer(serializers.ModelSerializer):

    class Meta:
        Model = Crystal
        fields = (
            'compound',
        )


