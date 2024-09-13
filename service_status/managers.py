from django.apps import apps
from django.db.models import F, Manager, QuerySet


class ServiceStateQueryset(QuerySet):
    def to_frontend(self):
        Service = apps.get_model("service_status", "Service")

        qs = Service.objects.annotate(
            id=F("service"),
            name=F("display_name"),
            state=F("last_state"),
        ).order_by("service")

        return qs


class ServiceStateDataManager(Manager):
    def get_queryset(self):
        return ServiceStateQueryset(self.model, using=self._db)

    def to_frontend(self):
        return self.get_queryset().to_frontend().values("id", "name", "state")
