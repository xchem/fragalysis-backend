from django.apps import apps
from django.db.models import F, Manager


class ServiceStateDataManager(Manager):
    def to_frontend(self):
        Service = apps.get_model("service_status", "Service")
        return Service.objects.order_by("service").values(
            id=F("service"),
            name=F("display_name"),
            state=F("last_state"),
        )
