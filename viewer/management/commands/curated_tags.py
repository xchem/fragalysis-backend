from django.core.management.base import BaseCommand

from viewer.utils import dump_curated_tags, restore_curated_tags


class Command(BaseCommand):
    help = "Dump or load curated tags"

    def add_arguments(self, parser):
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument(
            "--dump", metavar="<JSON file>", type=str, help="Save to file"
        )
        group.add_argument(
            "--load",
            metavar="JSON file>",
            type=str,
            help="Load file",
        )

    def handle(self, *args, **kwargs):
        # Unused args
        del args
        if kwargs["dump"]:
            dump_curated_tags(filename=kwargs["dump"])
        if kwargs["load"]:
            restore_curated_tags(filename=kwargs["load"])
