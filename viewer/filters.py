import django_filters
from django_filters import rest_framework as filters

from viewer.models import Snapshot
from viewer.models import SiteObservation
from viewer.models import CanonSite
from viewer.models import CanonSiteConf
from viewer.models import XtalformSite


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


class TargetFilterMixin(filters.FilterSet):
    target = filters.CharFilter(
        label="Target ID",
        field_name="target",
    )


class SiteObservationFilter(TargetFilterMixin):
    class Meta:
        model = SiteObservation
        fields = ("target",)


class CanonSiteFilter(TargetFilterMixin):
    class Meta:
        model = CanonSite
        fields = ("target",)


class CanonSiteConfFilter(TargetFilterMixin):
    class Meta:
        model = CanonSiteConf
        fields = ("target",)


class XtalformSiteFilter(TargetFilterMixin):
    class Meta:
        model = XtalformSite
        fields = ("target",)
