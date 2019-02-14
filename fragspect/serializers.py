from xchem_db.models import Crystal

from rest_framework import serializers



class FragspectSerializer(serializers.ModelSerializer):

    class Meta:
        Model = Crystal
        fields = (
            'compound',
        )


