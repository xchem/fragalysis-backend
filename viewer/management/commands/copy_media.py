"""Copy an existing stack's media to (or from) the S3 bucket (issue 1038).

The data migration from the filesystem (NFS) media storage to S3. Run it in a
backend Pod that has the media directory mounted and the MEDIA_S3_* settings.
It works whatever MEDIA_STORAGE is, so it can run while the stack still uses
the filesystem. See viewer/media_copy.py for what is (and is not) copied::

    python manage.py copy_media --dry-run           # counts and bytes, copies nothing
    python manage.py copy_media                     # copy (run again to resume)
    python manage.py copy_media --verify            # every file there, same size?
    python manage.py copy_media --check-references  # every DB file value there?
    python manage.py copy_media --direction s3-to-fs   # the other way (rollback)

A run copies only files that are missing (or have a different size) on the
target, so it can be repeated. --verify and --check-references copy nothing,
and fail (non-zero exit) if they find a problem.
"""
from pathlib import Path

import humanize
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured
from django.core.management.base import BaseCommand, CommandError

from frag_s3_storage import S3Storage
from viewer.media_copy import (
    FS_TO_S3,
    S3_TO_FS,
    copy_files,
    list_filesystem,
    list_s3,
    plan_copy,
    referenced_names,
)

# Longer listings are truncated
_MAX_LISTED = 20


class Command(BaseCommand):
    help = "Copy media between the media directory and the S3 bucket."

    def add_arguments(self, parser):
        parser.add_argument(
            "--direction",
            choices=(FS_TO_S3, S3_TO_FS),
            default=FS_TO_S3,
            help=f"The copy direction (default: {FS_TO_S3}).",
        )
        parser.add_argument(
            "--source",
            default=settings.MEDIA_ROOT,
            help="The media directory (default: MEDIA_ROOT).",
        )
        parser.add_argument(
            "--include-tmp", action="store_true", help="Also copy tmp/."
        )
        parser.add_argument(
            "--dry-run", action="store_true", help="Report only, copy nothing."
        )
        parser.add_argument(
            "--verify",
            action="store_true",
            help="Copy nothing. Fail if a file is missing or differs on the target.",
        )
        parser.add_argument(
            "--check-references",
            action="store_true",
            help="Copy nothing. Fail if a file that the database references is"
            " on the source but not on the target.",
        )

    def handle(self, *args, **kwargs):
        del args
        direction = kwargs["direction"]
        root = Path(kwargs["source"])
        try:
            storage = S3Storage(**settings.MEDIA_S3_STORAGE_OPTIONS)
        except ImproperlyConfigured as exc:
            raise CommandError(f"{exc} (set MEDIA_S3_BUCKET_NAME)") from exc

        filesystem = list_filesystem(root)
        bucket = list_s3(storage)
        source, target = (
            (filesystem, bucket) if direction == FS_TO_S3 else (bucket, filesystem)
        )
        plan = plan_copy(source, target, storage, direction, kwargs["include_tmp"])

        self.stdout.write(f"{direction}: {root} <-> s3://{storage.bucket_name}")
        for prefix, (count, size) in plan.by_prefix().items():
            self.stdout.write(
                f"  {prefix:<24} {count} files"
                f"  {humanize.naturalsize(size, binary=True)}"
            )
        self._print_list("skipped", plan.problems)

        if kwargs["check_references"]:
            self._check_references(source, target)
        elif kwargs["verify"]:
            self._verify(plan)
        elif kwargs["dry_run"]:
            size = sum(plan.wanted[name] for name in plan.to_copy)
            self.stdout.write(
                f"Dry run: would copy {len(plan.to_copy)} files"
                f" ({humanize.naturalsize(size, binary=True)})."
            )
        else:
            self._copy(plan, direction, root, storage)

    def _copy(self, plan, direction, root, storage) -> None:
        errors = copy_files(plan.to_copy, direction, root, storage)
        copied = len(plan.to_copy) - len(errors)
        self.stdout.write(self.style.SUCCESS(f"Done: copied {copied} files."))
        if errors:
            self._print_list("failed", errors)
            raise CommandError(f"{len(errors)} files were not copied. Run it again.")

    def _verify(self, plan) -> None:
        if plan.to_copy:
            self._print_list("missing or different on the target", plan.to_copy)
            raise CommandError(
                f"{len(plan.to_copy)} files are missing or different on the target."
            )
        self.stdout.write(self.style.SUCCESS(f"Verified {len(plan.wanted)} files."))

    def _check_references(self, source, target) -> None:
        referenced = referenced_names()
        not_on_target = referenced - set(target.files)
        not_on_source = referenced - set(source.files)
        lost = sorted(not_on_target - not_on_source)
        self.stdout.write(f"{len(referenced)} referenced files")
        self._print_list(
            "referenced, missing on the target (and missing on the source too)",
            sorted(not_on_target & not_on_source),
        )
        if lost:
            self._print_list("referenced, on the source but not on the target", lost)
            raise CommandError(
                f"{len(lost)} referenced file(s) are on the source but not on the"
                " target."
            )
        self.stdout.write(
            self.style.SUCCESS("Every referenced file on the source is on the target.")
        )

    def _print_list(self, title: str, lines: list[str]) -> None:
        if not lines:
            return
        self.stdout.write(f"  {title} ({len(lines)}):")
        for line in lines[:_MAX_LISTED]:
            self.stdout.write(f"    {line}")
        if len(lines) > _MAX_LISTED:
            self.stdout.write(f"    ... and {len(lines) - _MAX_LISTED} more")
