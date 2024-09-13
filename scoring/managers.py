import logging

from django.apps import apps
from django.db.models import F, Manager, QuerySet

logger = logging.getLogger(__name__)


class ScoreChoiceQueryset(QuerySet):
    def filter_qs(self):
        ScoreChoice = apps.get_model("scoring", "ScoreChoice")
        qs = ScoreChoice.objects.prefetch_related(
            "site_observation",
            "site_observation__experiment",
            "site_observation__experiment__experiment_upload",
            "site_observation__experiment__experiment_upload__target",
        ).annotate(
            target=F("site_observation__experiment__experiment_upload__target"),
        )

        return qs


class ScoreChoiceDataManager(Manager):
    def get_queryset(self):
        return ScoreChoiceQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class SiteObservationChoiceQueryset(QuerySet):
    def filter_qs(self):
        SiteObservationChoice = apps.get_model("scoring", "SiteObservationChoice")
        qs = SiteObservationChoice.objects.prefetch_related(
            "site_observation",
            "site_observation__experiment",
            "site_observation__experiment__experiment_upload",
            "site_observation__experiment__experiment_upload__target",
        ).annotate(
            target=F("site_observation__experiment__experiment_upload__target"),
        )

        return qs


class SiteObservationChoiceDataManager(Manager):
    def get_queryset(self):
        return SiteObservationChoiceQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class CmpdChoiceQueryset(QuerySet):
    def filter_qs(self):
        CmpdChoice = apps.get_model("scoring", "CmpdChoice")
        qs = CmpdChoice.objects.prefetch_related(
            "cmpd_id__siteobservation",
            "cmpd_id__siteobservation__experiment",
            "cmpd_id__siteobservation__experiment__experiment_upload",
            "cmpd_id__siteobservation__experiment__experiment_upload__target",
        ).annotate(
            target=F("cmpd_id__siteobservation__experiment__experiment_upload__target"),
        )

        return qs


class CmpdChoiceDataManager(Manager):
    def get_queryset(self):
        return CmpdChoiceQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class ViewSceneQueryset(QuerySet):
    def filter_qs(self):
        ViewScene = apps.get_model("scoring", "ViewScene")
        qs = ViewScene.objects.prefetch_related(
            "snapshot__session_project__target",
        ).annotate(
            target=F("snapshot__session_project__target"),
        )

        return qs


class ViewSceneDataManager(Manager):
    def get_queryset(self):
        return ViewSceneQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)


class SiteObservationAnnotationQueryset(QuerySet):
    def filter_qs(self):
        SiteObservationAnnotation = apps.get_model(
            "scoring", "SiteObservationAnnotation"
        )
        qs = SiteObservationAnnotation.objects.prefetch_related(
            "site_observation",
            "site_observation__experiment",
            "site_observation__experiment__experiment_upload",
            "site_observation__experiment__experiment_upload__target",
        ).annotate(
            target=F("site_observation__experiment__experiment_upload__target"),
        )

        return qs


class SiteObservationAnnotationDataManager(Manager):
    def get_queryset(self):
        return SiteObservationChoiceQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()

    def by_target(self, target):
        return self.get_queryset().filter_qs().filter(target=target.id)
