"""Copy media between a filesystem directory and the S3 bucket (issue 1038).

This is the core of the ``copy_media`` command: the data migration that moves an
existing stack's media from its (NFS) media directory to S3, or back. The S3 side
is always the bucket of MEDIA_S3_STORAGE_OPTIONS, whatever MEDIA_STORAGE is, so
the command can run on a stack that still uses the filesystem.

Object keys are the files' paths relative to the media directory, so the layout
does not change and a rollback is a copy in the other direction. Each side is
listed once, and only files that are missing or have a different size are
copied (in parallel), so a run can be repeated to resume, or to copy the
changes at the cutover.

Not copied: tmp/ (unless asked), squonk_upload/ (the Squonk2 integration is
deprecated), the media deletion watcher's files, and download archives that can
no longer be fetched. Symlinks, names that are not UTF-8 and names that are too
long for an S3 key are reported and skipped.
"""
import os
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

from django.apps import apps
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import SuspiciousFileOperation
from django.db import models

from frag_s3_storage import S3Storage
from viewer import media_storage
from viewer.media_watcher import PAUSE_PREFIX
from viewer.models import DownloadLinks

FS_TO_S3 = "fs-to-s3"
S3_TO_FS = "s3-to-fs"

# The longest object key that S3 accepts (UTF-8 bytes)
_MAX_KEY_BYTES = 1024

_DELETIONS_LOG = re.compile(
    rf"^{re.escape(settings.TARGET_LOADER_MEDIA_DIRECTORY)}/deletions-[^/]*\.log$"
)


@dataclass
class Listing:
    """The files on one side: {name: size}, and the entries that were skipped."""

    files: dict[str, int] = field(default_factory=dict)
    problems: list[str] = field(default_factory=list)


@dataclass
class CopyPlan:
    """What to copy from the source side to the target side."""

    # The source files that are wanted, as {name: size}
    wanted: dict[str, int] = field(default_factory=dict)
    # The wanted files that are missing (or have a different size) on the target
    to_copy: list[str] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    def by_prefix(self) -> dict[str, tuple[int, int]]:
        """{top-level prefix: (file count, bytes)} of the wanted files."""
        totals: dict[str, list[int]] = defaultdict(lambda: [0, 0])
        for name, size in self.wanted.items():
            total = totals[name.split("/", 1)[0]]
            total[0] += 1
            total[1] += size
        return {prefix: (n, size) for prefix, (n, size) in sorted(totals.items())}


def list_filesystem(root: Path) -> Listing:
    """Lists the files below a local directory (symlinks are not followed)."""
    listing = Listing()
    for directory, dir_names, file_names in os.walk(root, followlinks=False):
        for dir_name in dir_names:
            if Path(directory, dir_name).is_symlink():
                listing.problems.append(f"{Path(directory, dir_name)}: a symlink")
        for file_name in file_names:
            path = Path(directory, file_name)
            if path.is_symlink():
                listing.problems.append(f"{path}: a symlink")
                continue
            name = path.relative_to(root).as_posix()
            try:
                name.encode("utf-8")
            except UnicodeEncodeError:
                listing.problems.append(f"{os.fsencode(path)!r}: the name is not UTF-8")
                continue
            listing.files[name] = path.lstat().st_size
    return listing


def list_s3(storage: S3Storage) -> Listing:
    """Lists the objects in the bucket."""
    return Listing(files={name: size for name, size, _ in storage.iter_objects()})


def fetchable_downloads() -> set[str]:
    """The download archives that can still be fetched: those of static links,
    and of dynamic links that have not expired (and are not erased)."""
    links = (
        DownloadLinks.objects.filter(deleted__isnull=True)
        .exclude(file_url__isnull=True)
        .exclude(task_id__isnull=True)
        .filter(models.Q(static_link=True) | models.Q(expired_date__isnull=True))
    )
    return {name for link in links if (name := link.get_file_url())}


def _wanted(name: str, include_tmp: bool, downloads: set[str]) -> bool:
    """True if a source file is to be copied."""
    top = name.split("/", 1)[0]
    if top == "tmp":
        return include_tmp
    if top == "squonk_upload":
        return False
    if top == "downloads":
        return name in downloads
    if PurePosixPath(name).name.startswith(PAUSE_PREFIX):
        return False
    return not _DELETIONS_LOG.match(name)


def _s3_key(storage: S3Storage, name: str) -> str:
    return f"{storage.location}/{name}" if storage.location else name


def plan_copy(
    source: Listing,
    target: Listing,
    storage: S3Storage,
    direction: str,
    include_tmp: bool,
) -> CopyPlan:
    """Works out which source files to copy to the target."""
    plan = CopyPlan(problems=list(source.problems))
    downloads = fetchable_downloads()
    for name, size in sorted(source.files.items()):
        if not _wanted(name, include_tmp, downloads):
            continue
        if (
            direction == FS_TO_S3
            and len(_s3_key(storage, name).encode("utf-8")) > _MAX_KEY_BYTES
        ):
            plan.problems.append(f"{name}: the S3 key is longer than 1024 bytes")
            continue
        plan.wanted[name] = size
        if target.files.get(name) != size:
            plan.to_copy.append(name)
    return plan


def copy_files(
    names: list[str], direction: str, root: Path, storage: S3Storage
) -> list[str]:
    """Copies the files, in parallel. Returns the errors (one for each failure)."""
    _ = storage.client  # Make the (thread-safe) client before the threads

    def _copy(name: str) -> str | None:
        local = root.joinpath(name)
        try:
            if direction == FS_TO_S3:
                storage.upload_file(local, name)
            else:
                local.parent.mkdir(parents=True, exist_ok=True)
                storage.download_file(name, local)
        except Exception as exc:  # pylint: disable=broad-except
            # One failed file must not stop the copy: it is reported, and the
            # next run copies it again.
            return f"{name}: {exc}"
        return None

    threads = max(1, settings.MEDIA_TRANSFER_THREADS)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        return [error for error in executor.map(_copy, names) if error]


def _referencing_fields(model) -> list[tuple[str, bool]]:
    """(field name, is a list) of the model's fields that hold media file names."""
    found = []
    for model_field in model._meta.get_fields():  # pylint: disable=protected-access
        if isinstance(model_field, ArrayField) and isinstance(
            model_field.base_field, models.FileField
        ):
            found.append((model_field.name, True))
        elif isinstance(model_field, models.FileField):
            found.append((model_field.name, False))
    if model.__name__ in ("ComputedSet", "HistoricalComputedSet"):
        found.append(("written_sdf_filename", False))  # A TextField
    return found


def referenced_names() -> set[str]:
    """The storage names of every media file that a row (or history row) holds,
    and of the download archives that can be fetched."""
    names: set[str] = set()
    for model in apps.get_models():
        for field_name, is_list in _referencing_fields(model):
            for value in model.objects.values_list(field_name, flat=True).iterator():
                for item in (value or []) if is_list else [value]:
                    text = str(item or "").strip()
                    if text in ("", "None"):
                        continue
                    try:
                        names.add(media_storage.to_name(text))
                    except SuspiciousFileOperation:
                        continue
    return names | fetchable_downloads()
