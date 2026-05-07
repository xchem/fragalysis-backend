from django.core.management.base import BaseCommand

from viewer.cache import clear_all_view_caches
from viewer.utils import change_target_project


class Command(BaseCommand):
    help = "Move target to a different project/proposal"

    def add_arguments(self, parser):
        parser.add_argument(
            "--target",
            type=str,
            help="Target's title (not display name!)",
            required=True,
        )
        parser.add_argument(
            "--current_project",
            type=str,
            help="Target's current project",
            required=True,
        )
        parser.add_argument(
            "--new_project",
            type=str,
            help="Move target to this project. Must already exist",
            required=True,
        )

    def handle(self, *args, **kwargs):
        # Unused args
        del args
        try:
            change_target_project(
                target_name=kwargs["target"],
                old_project_name=kwargs["current_project"],
                new_project_name=kwargs["new_project"],
            )
        except ValueError as exc:
            self.stdout.write(self.style.ERROR(exc.args[0]))
            return
        # Reproject changes which users see this target's data via
        # ISPyBSafeQuerySet's project filter; flush cached responses.
        clear_all_view_caches()
