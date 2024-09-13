from django.core.management.base import BaseCommand

from service_status.utils import init_services


class Command(BaseCommand):
    help = "Activate service health check queries on startup"

    def handle(self, *args, **kwargs):
        # Unused args
        del args, kwargs
        init_services()
