from xchem_db.models import Crystal
from xchem_db.serializers import FragspectSerializer

from api.security import ISpyBSafeQuerySet


class CrystalView(ISpyBSafeQuerySet):
    queryset = Crystal.objects.filter()
    filter_permissions = "visit__proposal"
    serializer_class = FragspectSerializer
    filter_fields = ("target__target_name",)

