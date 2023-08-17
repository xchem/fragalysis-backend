import django_filters
from django_filters import rest_framework as filters

from viewer.models import Snapshot
from viewer.models import SiteObservation


class SnapshotFilter(filters.FilterSet):
    session_project = django_filters.CharFilter(
        field_name="session_project", lookup_expr="id"
    )
    session_project__isnull = django_filters.BooleanFilter(
        field_name="session_project", lookup_expr="isnull"
    )

    class Meta:
        model = Snapshot
        fields = [
            "id",
            "type",
            "title",
            "author",
            "description",
            "created",
            "data",
            "session_project",
            "parent",
            "children",
        ]


class SiteObservationFilter(filters.FilterSet):
    target_name = filters.CharFilter(
        label="Target",
        field_name="target_name",
    )

    class Meta:
        model = SiteObservation
        fields = ("target_name",)
