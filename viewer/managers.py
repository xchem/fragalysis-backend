from django.apps import apps

from django.db.models import QuerySet
from django.db.models import Manager
from django.db.models import F


class SiteObservationQueryset(QuerySet):
    def filter_qs(self):
        SiteObservation = apps.get_model("viewer", "SiteObservation")
        qs = SiteObservation.objects.prefetch_related(
            "experiment",
            "experiment__experiment_upload",
            "experiment__experiment_upload__target",
        ).annotate(
            target=F("experiment__experiment_upload__target"),
            target_name=F("experiment__experiment_upload__target__title"),
        )

        return qs


class SiteObservationDataManager(Manager):
    def get_queryset(self):
        return SiteObservationQueryset(self.model, using=self._db)

    def filter_qs(self):
        return self.get_queryset().filter_qs()
