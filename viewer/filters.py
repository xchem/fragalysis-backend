import logging

import django_filters
from django.core.exceptions import ValidationError
from django_filters import rest_framework as filters
from pgvector.django import L2Distance

from viewer.models import (
    AtomCoordinates,
    CanonSite,
    CanonSiteConf,
    Compound,
    ComputedSet,
    Experiment,
    ExperimentUpload,
    Pose,
    QuatAssembly,
    Result,
    SiteObservation,
    SiteObservationQualityStatus,
    Snapshot,
    SnapshotScreenshot,
    Target,
    XtalformSite,
)

logger = logging.getLogger(__name__)


class SnapshotFilter(filters.FilterSet):
    session_project = django_filters.CharFilter(
        field_name="session_project", lookup_expr="id"
    )
    session_project__isnull = django_filters.BooleanFilter(
        field_name="session_project", lookup_expr="isnull"
    )

    target = django_filters.CharFilter(
        field_name="session_project__target",
        lookup_expr="id",
        label="Target",
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


class SnapshotScreenshotFilter(filters.FilterSet):
    class Meta:
        model = SnapshotScreenshot
        fields = [
            "snapshot",
            "screenshot_type",
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


class SiteObservationCoordinateFilter(TargetFilterMixin):
    # adding these as fields, but they're not in the model, hence the no-op method
    xorigin = django_filters.NumberFilter(method="noop", label="X origin")
    yorigin = django_filters.NumberFilter(method="noop", label="Y origin")
    zorigin = django_filters.NumberFilter(method="noop", label="Z origin")
    radius = django_filters.NumberFilter(method="noop", label="Radius")

    class Meta:
        model = SiteObservation
        fields = ("target", "xorigin", "yorigin", "zorigin", "radius")

    # I feel like this isn't the correct method in FilterSet class to
    # override for this functionality. but given the sparse
    # documentation of django_filter I can't find a better one either

    # Update: ok this is clearly not the correct place to update. I
    # want to add coords to be the last step, so it would already be
    # filtered by target, otherwise cannot guarantee consistent
    # behaviour
    def filter_queryset(self, queryset):
        qs = super().filter_queryset(queryset)
        qs = self.filter_by_radius(queryset)
        return qs

    def noop(self, queryset, name, value):
        del name, value
        return queryset

    def filter_by_radius(self, queryset):
        params = self.data
        logger.debug('params: %s', params)
        logger.debug('queryset count: %s', queryset.count())

        x_str = params.get("xorigin", None)
        y_str = params.get("yorigin", None)
        z_str = params.get("zorigin", None)
        r_str = params.get("radius", None)
        target_id = params.get("target", None)

        if not all([x_str, y_str, z_str, r_str, target_id]):
            # none of the coordinate filter parameters defined, it's
            # clear the user is not trying to filter by coordinate
            return queryset

        try:
            x = float(x_str)
        except TypeError as exc:
            raise ValidationError("xorigin is not a valid float value") from exc

        try:
            y = float(y_str)
        except TypeError as exc:
            raise ValidationError("yorigin is not a valid float value") from exc

        try:
            z = float(z_str)
        except TypeError as exc:
            raise ValidationError("zorigin is not a valid float value") from exc

        try:
            r = float(r_str)
        except TypeError as exc:
            raise ValidationError("Radius is not a valid float value") from exc

        if r < 0:
            raise ValidationError("Radius must be non-negative")

        try:
            target = Target.objects.get(pk=target_id)
        except Target.DoesNotExist as exc:
            raise ValidationError(f"Target with pk {target_id} not found") from exc

        logger.debug('x: %s', x)
        logger.debug('y: %s', y)
        logger.debug('z: %s', z)
        logger.debug('r: %s', r)
        logger.debug('target: %s', target)

        # all params present and valid, continue to filter. when target is defined

        # fmt: off
        qs = SiteObservation.filter_manager.by_target(target).filter(
            pk__in=AtomCoordinates.objects.alias(
                distance=L2Distance('coords', [x, y, z]),
            ).filter(
                distance__lte=r,
            ).values(
                'site_observation',
            ),
        )
        # fmt: on

        logger.debug('filtered qs: %s', qs.count())

        return qs


class CanonSiteFilter(TargetFilterMixin):
    class Meta:
        model = CanonSite
        fields = ("target",)


class ExperimentFilter(TargetFilterMixin):
    class Meta:
        model = Experiment
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


class ComputedSetFilter(filters.FilterSet):
    project = django_filters.CharFilter(
        field_name="project",
        lookup_expr="icontains",
        label="Project",
    )

    class Meta:
        model = ComputedSet
        fields = ("name", "target", "project")


class ExperimentUploadFilter(filters.FilterSet):
    class Meta:
        model = ExperimentUpload
        fields = (
            "target",
            "project",
            "committer",
            "data_version_major",
            "data_version_minor",
        )


class SiteObservationQualityStatusFilter(filters.FilterSet):
    class Meta:
        model = SiteObservationQualityStatus
        fields = (
            "site_observation",
            "status",
            "user",
            "timestamp",
            "auto_assigned",
            "main_status",
        )


class ActivityResultFilter(filters.FilterSet):
    class Meta:
        model = Result
        fields = (
            "result_upload__target",
            "result_property__data_type__data_type",
            "result_property__result_property",
            "result_property__unit",
            "compound",
            "site_observation",
        )
