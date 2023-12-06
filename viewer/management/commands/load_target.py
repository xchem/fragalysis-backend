from django.core.management.base import BaseCommand

from django.db import IntegrityError

from viewer.target_loader import load_target


class Command(BaseCommand):
    help = 'Load target directly from media/ directory.'

    def add_arguments(self, parser):
        parser.add_argument('data_bundle', type=str, help='Data archive to be loaded')

    def handle(self, *args, **kwargs):
        # Unused args
        del args

        try:
            load_target(kwargs["data_bundle"])
            self.stdout.write(self.style.SUCCESS("Data imported"))
        except KeyError as err:
            self.stdout.write(self.style.ERROR(err.args[0]))
        except IntegrityError as err:
            self.stdout.write(self.style.ERROR(err.args[0]))
        except FileNotFoundError as err:
            self.stdout.write(self.style.ERROR(err.args[0]))
