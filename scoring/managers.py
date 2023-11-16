import logging

from django.apps import apps

from django.db.models import QuerySet
from django.db.models import Manager
from django.db.models import F


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
