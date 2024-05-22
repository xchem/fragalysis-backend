from django.core.management.base import BaseCommand

from service_status.utils import services


class Command(BaseCommand):
    help = "Activate/deactivate service health check queries defined in service_status/services.py"

    def add_arguments(self, parser):
        parser.add_argument(
            "--enable", metavar="Service IDs", nargs="*", help="Enable service queries"
        )
        parser.add_argument(
            "--disable",
            metavar="Service IDs",
            nargs="*",
            help="Disable service queries",
        )

    def handle(self, *args, **kwargs):
        # Unused args
        del args

        if "enable" not in kwargs.keys() and "disable" not in kwargs.keys():
            self.stdout.write(
                self.style.ERROR("One of '--enable' or '--disable' must be defined'")
            )
            return

        services(enable=kwargs["enable"], disable=kwargs["disable"])
