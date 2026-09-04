"""Report - and with ``--restore``, rebuild - ``DownloadLinks`` rows for download
archives orphaned under ``media/downloads/``.

Migration ``0153_alter_downloadlinks_proteins`` empties the table before altering
the ``proteins`` column. It was written in January 2026 on a long-lived branch and
only merged to staging (2026-07-08) and production (2026-07-23), so on those
instances every download record made before that deploy was wiped while the
archives stayed on disk. A static link handed out earlier now fails with
"file_url should only be provided for static files".

The original POST body is embedded in the README.md inside each archive and the
directory name is the task ID, so the rows can be rebuilt - see
``viewer.download_restore`` for exactly what is (and is not) restored.

Nothing is written without ``--restore``. Run it without arguments first and read
the report::

    python manage.py restore_download_links
    python manage.py restore_download_links --restore

Static links are restored by default. ``--include-dynamic`` also restores dynamic
links with their original (long expired) retention, which hands them to the normal
housekeeping and so **deletes those archives** after the hard-expiry grace period.
``--unparsed-as-static`` covers archives built with ``readme=False``, whose search
is unrecoverable: the link is restored so the URL works, with no search recorded.
"""

import humanize
from django.conf import settings
from django.core.management.base import BaseCommand

from viewer.download_restore import RestoreReport, restore, scan_downloads
from viewer.models import Target

# Longer listings are truncated unless --list-all is given.
_MAX_LISTED = 100


class Command(BaseCommand):
    help = "Report orphaned download archives; --restore rebuilds their DB records."

    def add_arguments(self, parser):
        parser.add_argument(
            "--restore",
            action="store_true",
            help="Actually create the records. Without this nothing is written.",
        )
        parser.add_argument(
            "--include-dynamic",
            action="store_true",
            help="Also restore dynamic (non-static) links. Their retention has long"
            " passed, so the housekeeping will expire them and delete the archives.",
        )
        parser.add_argument(
            "--unparsed-as-static",
            action="store_true",
            help="Also restore archives with no readable download command, as static"
            " links with no search recorded.",
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
        do_restore = kwargs.get("restore", False)
        include_dynamic = kwargs.get("include_dynamic", False)
        unparsed_as_static = kwargs.get("unparsed_as_static", False)
        force = kwargs.get("force", False)
        list_all = kwargs.get("list_all", False)

        if not settings.MEDIA_ROOT:
            self.stdout.write(self.style.ERROR("MEDIA_ROOT is not set."))
            return

        report = scan_downloads()
        self._print_report(report, list_all=list_all)

        selected = self._selected(report, include_dynamic, unparsed_as_static)

        if not do_restore:
            self.stdout.write("")
            self.stdout.write(
                self.style.SUCCESS(
                    "Report only - nothing written. Re-run with --restore to"
                    f" rebuild {len(selected)} records."
                )
            )
            return

        # Everything on disk looking orphaned is the signature of running against
        # the wrong database - the same guard cleanup_media applies.
        if not Target.objects.exists() and not force:
            self.stdout.write(
                self.style.ERROR(
                    "The database holds no targets at all. Refusing to write"
                    " - re-run with --force if that is really the case."
                )
            )
            return

        if include_dynamic:
            self.stdout.write(
                self.style.WARNING(
                    "--include-dynamic: restored dynamic links are already past"
                    " their retention, so the download cleanup will expire them"
                    " and delete their archives."
                )
            )

        result = restore(
            report.orphans,
            include_dynamic=include_dynamic,
            unparsed_as_static=unparsed_as_static,
        )

        self.stdout.write("")
        self.stdout.write(
            self.style.SUCCESS(
                f"Restored {len(result.restored)} records"
                f" ({humanize.naturalsize(result.restored_bytes, binary=True)}"
                " of archives back in service)."
            )
        )
        if result.failed:
            self.stdout.write(self.style.ERROR(f"  {len(result.failed)} failed:"))
            for orphan, message in result.failed[:_MAX_LISTED]:
                self.stdout.write(f"    {orphan.task_id}  {message}")

    @staticmethod
    def _selected(report: RestoreReport, include_dynamic: bool, unparsed: bool) -> list:
        selected = list(report.static_orphans)
        if include_dynamic:
            selected += report.dynamic_orphans
        if unparsed:
            selected += report.unparsed_orphans
        return selected

    def _print_report(self, report: RestoreReport, list_all: bool) -> None:
        self.stdout.write("")
        self.stdout.write(f"downloads  ({report.root})")
        self.stdout.write(f"  scanned                {report.scanned}")
        self.stdout.write(f"  referenced             {report.referenced}")
        self.stdout.write(f"  orphaned               {len(report.orphans)}")
        self.stdout.write(f"    static               {len(report.static_orphans)}")
        self.stdout.write(f"    dynamic              {len(report.dynamic_orphans)}")
        self.stdout.write(f"    no search recorded   {len(report.unparsed_orphans)}")
        self.stdout.write(f"  skipped                {len(report.skipped)}")

        self._print_list(
            "static (restored by default)",
            [self._describe(o) for o in report.static_orphans],
            list_all,
        )
        self._print_list(
            "dynamic (--include-dynamic)",
            [self._describe(o) for o in report.dynamic_orphans],
            list_all,
        )
        self._print_list(
            "no search recorded (--unparsed-as-static)",
            [self._describe(o) for o in report.unparsed_orphans],
            list_all,
        )
        self._print_list(
            "skipped entries",
            [f"{path.name}  ({reason})" for path, reason in report.skipped],
            list_all,
        )
        self._print_list("warnings", report.warnings, list_all)

    @staticmethod
    def _describe(orphan) -> str:
        target = orphan.target.title if orphan.target else "?"
        return (
            f"{orphan.task_id}/{orphan.file_url}"
            f"  {humanize.naturalsize(orphan.size_bytes, binary=True)}"
            f"  {orphan.create_date:%Y-%m-%d}  target={target}"
        )

    def _print_list(self, title: str, lines: list[str], list_all: bool) -> None:
        if not lines:
            return
        self.stdout.write(f"  {title}:")
        shown = lines if list_all else lines[:_MAX_LISTED]
        for line in shown:
            self.stdout.write(f"    {line}")
        if len(lines) > len(shown):
            self.stdout.write(f"    ... and {len(lines) - len(shown)} more")
