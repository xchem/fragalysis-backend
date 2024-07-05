import logging

from django.apps import apps
from django.db.models import F, Manager, OuterRef, QuerySet, Subquery

logger = logging.getLogger(__name__)


class SiteObservationQueryset(QuerySet):
    def filter_qs(self):
        SiteObservation = apps.get_model("viewer", "SiteObservation")
        qs = SiteObservation.objects.prefetch_related(
            "experiment",
            "experiment__experiment_upload",
            "experiment__experiment_upload__target",
            "cmpd",
        ).annotate(
            target=F("experiment__experiment_upload__target"),
            compound_code=F("cmpd__compound_code"),
            prefix_tooltip=F("experiment__prefix_tooltip"),
        )

        return qs


class SiteObservationDataManager(Manager):
    def get_queryset(self):
        return SiteObservationQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class ExperimentQueryset(QuerySet):
    def filter_qs(self):
        Experiment = apps.get_model("viewer", "Experiment")
        qs = Experiment.objects.prefetch_related(
            "experiment_upload",
            "experiment_upload__target",
        ).annotate(
            target=F("experiment_upload__target"),
        )

        return qs


class ExperimentDataManager(Manager):
    def get_queryset(self):
        return ExperimentQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class CompoundQueryset(QuerySet):
    def filter_qs(self):
        Compound = apps.get_model("viewer", "Compound")
        SiteObservation = apps.get_model("viewer", "SiteObservation")

        so_qs = SiteObservation.filter_manager.filter_qs()

        qs = Compound.objects.annotate(
            target=Subquery(
                so_qs.filter(
                    cmpd=OuterRef("pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class CompoundDataManager(Manager):
    def get_queryset(self):
        return CompoundQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.pk)


class XtalformQueryset(QuerySet):
    def filter_qs(self):
        Xtalform = apps.get_model("viewer", "Xtalform")
        Experiment = apps.get_model("viewer", "Experiment")

        exp_qs = Experiment.filter_manager.filter_qs()

        qs = Xtalform.objects.annotate(
            target=Subquery(
                exp_qs.filter(
                    xtalform=OuterRef("pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class XtalformDataManager(Manager):
    def get_queryset(self):
        return XtalformQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class XtalformSiteQueryset(QuerySet):
    def filter_qs(self):
        XtalformSite = apps.get_model("viewer", "XtalformSite")
        Experiment = apps.get_model("viewer", "Experiment")

        exp_qs = Experiment.filter_manager.filter_qs()

        qs = XtalformSite.objects.prefetch_related(
            "xtalform",
        ).annotate(
            target=Subquery(
                exp_qs.filter(
                    xtalform=OuterRef("xtalform__pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class XtalformSiteDataManager(Manager):
    def get_queryset(self):
        return XtalformSiteQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


# this is perhaps a bit unusual, a dedicated mangager for m2m junction
# table. but the contents are coming from yaml file and are used to
# determin related objects (QuatAssembly) so I think in this case this
# is justified
class XtalformQuatAssemblyQueryset(QuerySet):
    def filter_qs(self):
        Xtalform = apps.get_model("viewer", "Xtalform")
        XtalformQuatAssembly = apps.get_model("viewer", "XtalformQuatAssembly")

        xtal_qs = Xtalform.filter_manager.filter_qs()

        qs = XtalformQuatAssembly.objects.prefetch_related(
            "xtalform",
        ).annotate(
            target=Subquery(
                xtal_qs.filter(
                    pk=OuterRef("xtalform__pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class XtalformQuatAssemblyDataManager(Manager):
    def get_queryset(self):
        return XtalformQuatAssemblyQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class QuatAssemblyQueryset(QuerySet):
    def filter_qs(self):
        QuatAssembly = apps.get_model("viewer", "QuatAssembly")
        XtalformQuatAssembly = apps.get_model("viewer", "XtalformQuatAssembly")

        xtlq_qs = XtalformQuatAssembly.filter_manager.filter_qs()

        qs = QuatAssembly.objects.annotate(
            target=Subquery(
                xtlq_qs.filter(
                    quat_assembly=OuterRef("pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class QuatAssemblyDataManager(Manager):
    def get_queryset(self):
        return QuatAssemblyQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class CanonSiteQueryset(QuerySet):
    def filter_qs(self):
        XtalformSite = apps.get_model("viewer", "XtalformSite")
        CanonSite = apps.get_model("viewer", "CanonSite")

        xtalsite_qs = XtalformSite.filter_manager.filter_qs()

        qs = CanonSite.objects.annotate(
            target=Subquery(
                xtalsite_qs.filter(
                    canon_site=OuterRef("pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class CanonSiteDataManager(Manager):
    def get_queryset(self):
        return CanonSiteQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class CanonSiteConfQueryset(QuerySet):
    def filter_qs(self):
        # this could equally be filtered via siteobservation.. any
        # advantage in doing that?
        CanonSite = apps.get_model("viewer", "CanonSite")
        CanonSiteConf = apps.get_model("viewer", "CanonSiteConf")

        canonsite_qs = CanonSite.filter_manager.filter_qs()

        qs = CanonSiteConf.objects.annotate(
            target=Subquery(
                canonsite_qs.filter(
                    pk=OuterRef("canon_site__pk"),
                ).values(
                    "target"
                )[:1],
            ),
        )

        return qs


class CanonSiteConfDataManager(Manager):
    def get_queryset(self):
        return CanonSiteConfQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class PoseQueryset(QuerySet):
    def filter_qs(self):
        Pose = apps.get_model("viewer", "Pose")
        CanonSite = apps.get_model("viewer", "CanonSite")

        canonsite_qs = CanonSite.filter_manager.filter_qs()

        qs = Pose.objects.annotate(
            target=Subquery(
                canonsite_qs.filter(
                    pk=OuterRef("canon_site__pk"),
                ).values(
                    "target"
                )[:1],
            ),
            upload_name=F("main_site_observation__code"),
        )

        return qs


class PoseDataManager(Manager):
    def get_queryset(self):
        return PoseQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class SnapshotQueryset(QuerySet):
    def filter_qs(self):
        Snapshot = apps.get_model("viewer", "Snapshot")
        qs = Snapshot.objects.prefetch_related(
            "session_project",
            "session_project__target",
        ).annotate(
            target=F("session_project__target"),
        )

        return qs


class SnapshotDataManager(Manager):
    def get_queryset(self):
        return SnapshotQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class SessionActionsQueryset(QuerySet):
    def filter_qs(self):
        SessionActions = apps.get_model("viewer", "SessionActions")
        qs = SessionActions.objects.prefetch_related(
            "session_project",
            "session_project__target",
        ).annotate(
            target=F("session_project__target"),
        )

        return qs


class SessionActionsDataManager(Manager):
    def get_queryset(self):
        return SessionActionsQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class CompoundIdentifierQueryset(QuerySet):
    def filter_qs(self):
        CompoundIdentifier = apps.get_model("viewer", "CompoundIdentifier")
        qs = CompoundIdentifier.objects.prefetch_related(
            "cmpd_id__siteobservation",
            "cmpd_id__siteobservation__experiment",
            "cmpd_id__siteobservation__experiment__experiment_upload",
            "cmpd_id__siteobservation__experiment__experiment_upload__target",
        ).annotate(
            target=F("cmpd_id__siteobservation__experiment__experiment_upload__target"),
        )

        return qs


class CompoundIdentifierDataManager(Manager):
    def get_queryset(self):
        return CompoundIdentifierQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)
