"""Report - and with ``--delete``, remove - orphaned files in the media directory.

Two areas are covered, each cross-referenced against the database by
``viewer.media_cleanup``:

* ``media/target_loader_data/`` - directories no ``Target.zip_archive`` names.
  The bulk of the debris: before target deletion was a feature, targets were
  removed by running commands in the Django shell, which left the whole
  per-target directory behind.
* ``media/computed_set_data/`` - files no computed-set row points at.

Nothing is removed without ``--delete``. Run it without arguments first and read
the report::

    python manage.py cleanup_media
    python manage.py cleanup_media --delete

Symlinks, non-directories in ``target_loader_data/`` (the file watcher writes its
log there) and anything resolving outside the area are never touched, and every
removal is logged.
"""

import humanize
from django.conf import settings
from django.core.management.base import BaseCommand

from viewer.media_cleanup import (
    AreaReport,
    delete_candidates,
    scan_computed_set_data,
    scan_target_loader_data,
    unexpected_computed_set_fields,
)
from viewer.models import ComputedSet, SiteObservation, Target

# Longer listings are truncated unless --list-all is given.
_MAX_LISTED = 100


class Command(BaseCommand):
    help = "Report orphaned media files and directories; --delete removes them."

    def add_arguments(self, parser):
        parser.add_argument(
            "--delete",
            action="store_true",
            help="Actually remove the orphans. Without this nothing is written.",
        )
        parser.add_argument(
            "--target-loader-data",
            action="store_true",
            help="Only scan target_loader_data/ (default: both areas).",
        )
        parser.add_argument(
            "--computed-set-data",
            action="store_true",
            help="Only scan computed_set_data/ (default: both areas).",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Skip the empty-database guard. Nothing else is relaxed.",
        )
        parser.add_argument(
            "--list-all",
            action="store_true",
            help=f"List every entry instead of the first {_MAX_LISTED}.",
        )

    def handle(self, *args, **kwargs):
        del args
        delete = kwargs.get("delete", False)
        force = kwargs.get("force", False)
        list_all = kwargs.get("list_all", False)

        only_loader = kwargs.get("target_loader_data", False)
        only_computed = kwargs.get("computed_set_data", False)
        # Neither flag means both areas.
        do_loader = only_loader or not only_computed
        do_computed = only_computed or not only_loader

        media_root = settings.MEDIA_ROOT
        if not media_root:
            self.stdout.write(self.style.ERROR("MEDIA_ROOT is not set."))
            return

        # A new FileField pointing at computed_set_data/ that media_cleanup
        # doesn't know about would make its files look orphaned.
        unexpected = unexpected_computed_set_fields()
        if unexpected:
            self.stdout.write(
                self.style.ERROR(
                    "Unknown fields point at computed_set_data/: "
                    + ", ".join(unexpected)
                    + ". Add them to _COMPUTED_SET_FIELDS in viewer/media_cleanup.py"
                    " before deleting anything."
                )
            )

        reports = []
        if do_loader:
            reports.append(scan_target_loader_data())
        if do_computed:
            reports.append(scan_computed_set_data())

        for report in reports:
            self._print_report(report, list_all=list_all)

        total_bytes = sum(r.reclaimable_bytes for r in reports)
        total_count = sum(len(r.candidates) for r in reports)

        if not delete:
            self.stdout.write("")
            self.stdout.write(
                self.style.SUCCESS(
                    f"Report only - nothing removed. Re-run with --delete to"
                    f" remove {total_count} entries"
                    f" ({humanize.naturalsize(total_bytes, binary=True)})."
                )
            )
            return

        if unexpected:
            self.stdout.write(
                self.style.ERROR("Refusing to delete while unknown fields exist.")
            )
            return

        for report in reports:
            if guard := self._empty_database_guard(report):
                if not force:
                    self.stdout.write(
                        self.style.ERROR(
                            f"{report.area}: {guard} Refusing to delete"
                            " - re-run with --force if that is really the case."
                        )
                    )
                    continue
                self.stdout.write(f"{report.area}: {guard} Proceeding (--force).")

            deleted, failed, errors = delete_candidates(report)
            self.stdout.write("")
            self.stdout.write(
                self.style.SUCCESS(
                    f"{report.area}: removed {deleted} entries"
                    f" ({humanize.naturalsize(report.reclaimable_bytes, binary=True)})."
                )
            )
            if failed:
                self.stdout.write(self.style.ERROR(f"  {failed} could not be removed:"))
                for message in errors[:_MAX_LISTED]:
                    self.stdout.write(f"    {message}")

    @staticmethod
    def _empty_database_guard(report: AreaReport) -> str | None:
        """A warning if the area's tables are empty, meaning everything on disk
        looks orphaned - the signature of running against the wrong database."""
        if report.area == settings.TARGET_LOADER_MEDIA_DIRECTORY:
            if not Target.objects.exists():
                return "the database holds no targets at all."
        elif not (ComputedSet.objects.exists() or SiteObservation.objects.exists()):
            return "the database holds no computed sets or site observations at all."
        return None

    def _print_report(self, report: AreaReport, list_all: bool) -> None:
        self.stdout.write("")
        self.stdout.write(f"{report.area}  ({report.root})")
        self.stdout.write(f"  scanned                {report.scanned}")
        self.stdout.write(f"  referenced             {report.referenced}")
        self.stdout.write(
            f"  orphaned               {len(report.candidates)}"
            f"   {humanize.naturalsize(report.reclaimable_bytes, binary=True)}"
        )
        self.stdout.write(f"  skipped                {len(report.skipped)}")
        self.stdout.write(f"  referenced but absent  {len(report.missing)}")

        self._print_list(
            "orphans",
            [
                f"{c.path.name}  {humanize.naturalsize(c.size_bytes, binary=True)}"
                for c in report.candidates
            ],
            list_all,
        )
        self._print_list(
            "skipped entries",
            [f"{path.name}  ({reason})" for path, reason in report.skipped],
            list_all,
        )
        self._print_list("referenced but absent", report.missing, list_all)
        self._print_list("warnings", report.warnings, list_all)

    def _print_list(self, title: str, lines: list[str], list_all: bool) -> None:
        if not lines:
            return
        self.stdout.write(f"  {title}:")
        shown = lines if list_all else lines[:_MAX_LISTED]
        for line in shown:
            self.stdout.write(f"    {line}")
        if len(lines) > len(shown):
            self.stdout.write(f"    ... and {len(lines) - len(shown)} more")
