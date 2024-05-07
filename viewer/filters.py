import django_filters
from django_filters import rest_framework as filters

from viewer.models import (
    CanonSite,
    CanonSiteConf,
    Compound,
    Pose,
    QuatAssembly,
    SiteObservation,
    Snapshot,
    XtalformSite,
)


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


class VectorFilter(TargetFilterMixin):
    class Meta:
        model = SiteObservation
        fields = ("id", "target", "cmpd", "smiles", "site_observation_groups")


class GraphFilter(TargetFilterMixin):
    class Meta:
        model = SiteObservation
        fields = ("target", "cmpd", "smiles", "site_observation_groups")


class MolpropsFilter(TargetFilterMixin):
    class Meta:
        model = Compound
        fields = ("target", "smiles", "inchi")


class MolImgFilter(TargetFilterMixin):
    class Meta:
        model = SiteObservation
        fields = ("target", "cmpd", "smiles", "site_observation_groups")


class CmpdImgFilter(TargetFilterMixin):
    class Meta:
        model = Compound
        fields = ("target", "smiles")


class CompoundFilter(TargetFilterMixin):
    class Meta:
        model = Compound
        fields = ("smiles", "current_identifier", "inchi")


class PoseFilter(TargetFilterMixin):
    class Meta:
        model = Pose
        fields = ("target", "canon_site", "compound", "main_site_observation")


class AssemblyFilter(TargetFilterMixin):
    class Meta:
        model = QuatAssembly
        fields = ("target",)
