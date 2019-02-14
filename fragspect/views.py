from xchem_db.models import Crystal
from xchem_db.serializers import CrystalSerializer

from api.security import ISpyBSafeQuerySet


class CrystalView(ISpyBSafeQuerySet):
    queryset = Crystal.objects.filter()
    filter_permissions = "visit__proposal"
    serializer_class = CrystalSerializer
    filter_fields = ("target__target_name",)

