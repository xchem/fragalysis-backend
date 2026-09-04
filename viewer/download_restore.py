"""Rebuild ``DownloadLinks`` rows for download archives orphaned in the media
directory.

Why this is needed
------------------
Migration ``0153_alter_downloadlinks_proteins`` empties the table
(``DownloadLinks.objects.all().delete()``) before changing the ``proteins``
column type. It was authored on 2026-01-12 on a long-lived feature branch and
only reached staging (``ccdc3ced``, 2026-07-08) and production (``7989f44e``,
2026-07-23), so on those instances the wipe ran at *that* deploy - taking out
every download record made up to then while leaving every archive under
``media/downloads/`` untouched on disk.

The visible symptom is a static download link handed out before that deploy
now failing with *"file_url should only be provided for static files"*:
``DownloadStructuresView.create()`` looks the URL up by ``task_id`` +
``file_url``, finds nothing, and falls through to that message.

Why the archives are recoverable
--------------------------------
``DownloadStructures._build_readme()`` writes the original POST body into the
README.md inside every archive (as a single fenced JSON line), and the
per-download directory name *is* the ``task_id``. Between them they carry
everything the lookup needs, so a row can be rebuilt from the archive alone.

What counts as an orphan
------------------------
``<MEDIA_ROOT>/downloads/<task_id>/<target>.{zip,tar.gz}``. Nothing in the code
deletes ``DownloadLinks`` rows other than migration 0153 and the ``Target`` /
``User`` cascades, so a directory whose name matches no ``DownloadLinks.task_id``
is a lost record. Directories that *are* referenced are left strictly alone, which
makes the restore idempotent.

What is restored, and what is deliberately not
----------------------------------------------
``task_id``, ``file_url``, ``static_link``, ``target``, ``create_date`` (the
archive's mtime - the request time is unknowable), ``original_search`` and the
``protein_params`` / ``other_params`` derived from it by the same
``get_download_params()`` the view uses.

``proteins`` is left NULL, on purpose. It holds ``SiteObservation`` *primary
keys* resolved at request time (see ``DownloadStructuresView.create()``); after
a re-upload those pks mean something else entirely. Leaving it NULL also
guarantees a restored row can never satisfy the reuse query in that view - which
matches on ``proteins=<a real string>`` - so no user can be silently handed a
stale archive. Restored rows serve the direct ``file_url`` lookup only.

``user`` is left NULL for the same reason: it isn't in the archive.

Dynamic links are excluded unless asked for. A dynamic record is restored with
its original retention (``create_date + KEEP_UNTIL_DURATION``), which is long
past - so the housekeeping in ``viewer.download_structures`` will expire it and
then delete its files, which is what should have happened months ago. That is a
real deletion, hence the separate opt-in.
"""

import json
import logging
import os
import re
import tarfile
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

from django.conf import settings
from django.db import transaction

from viewer.download_structures import KEEP_UNTIL_DURATION, get_download_params
from viewer.models import DownloadLinks, Target
from viewer.serializers import DownloadStructuresSerializer

logger = logging.getLogger(__name__)

# The subdirectory of MEDIA_ROOT holding one directory per download, named
# after the task that built it - see DownloadLinks.get_file_url().
DOWNLOADS_SUBDIR = "downloads"

# Archive suffixes create_download() can produce (use_zip switches between them).
_ARCHIVE_SUFFIXES = (".zip", ".tar.gz")

# The README line written by DownloadStructures._build_readme():
#   ```{"target_name": ..., ...}```
_README_JSON = re.compile(r"^```(\{.*\})```\s*$", re.MULTILINE)

# The README is a small text file; cap the read so a corrupt or hostile archive
# can't be expanded into memory.
_MAX_README_BYTES = 4 * 1024 * 1024


@dataclass(frozen=True)
class Orphan:
    """One download directory with no ``DownloadLinks`` row."""

    directory: Path
    archive: Path
    task_id: str
    file_url: str
    size_bytes: int
    create_date: datetime
    # None when the archive carries no readable original search (it was built
    # with readme=False, or the README could not be parsed).
    original_search: dict | None
    static_link: bool
    target: Target | None
    note: str | None = None

    @property
    def parsed(self) -> bool:
        return self.original_search is not None


@dataclass
class RestoreReport:
    """What a scan of ``media/downloads/`` found."""

    root: Path
    scanned: int = 0
    referenced: int = 0
    orphans: list[Orphan] = field(default_factory=list)
    skipped: list[tuple[Path, str]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def static_orphans(self) -> list[Orphan]:
        return [o for o in self.orphans if o.parsed and o.static_link]

    @property
    def dynamic_orphans(self) -> list[Orphan]:
        return [o for o in self.orphans if o.parsed and not o.static_link]

    @property
    def unparsed_orphans(self) -> list[Orphan]:
        return [o for o in self.orphans if not o.parsed]


@dataclass
class RestoreResult:
    restored: list[Orphan] = field(default_factory=list)
    failed: list[tuple[Orphan, str]] = field(default_factory=list)

    @property
    def restored_bytes(self) -> int:
        return sum(o.size_bytes for o in self.restored)


def downloads_root() -> Path:
    return Path(settings.MEDIA_ROOT).joinpath(DOWNLOADS_SUBDIR)


def _find_archive(directory: Path) -> tuple[Path | None, str | None]:
    """The single download archive in ``directory``, or a reason it can't be used."""
    try:
        entries = sorted(directory.iterdir())
    except OSError as exc:
        return None, f"unreadable ({exc})"

    archives = [
        entry
        for entry in entries
        if entry.is_file()
        and not entry.is_symlink()
        and entry.name.endswith(_ARCHIVE_SUFFIXES)
    ]
    if not archives:
        return None, "no .zip/.tar.gz archive"
    if len(archives) > 1:
        names = ", ".join(a.name for a in archives)
        return None, f"more than one archive ({names})"
    return archives[0], None


def _readme_text(archive: Path) -> str | None:
    """The README.md from the archive root, or None if there isn't one.

    Both archive flavours are written from ``.`` (``tar -C <dir> -cf - .`` and
    ``7z a -tzip <file> .``), so the member is ``./README.md`` or ``README.md``.
    """

    def _shallowest(names: list[str]) -> str | None:
        candidates = [
            name
            for name in names
            if PurePosixPath(name).name == "README.md"
            and len(PurePosixPath(name).parts) <= 2
        ]
        return min(candidates, key=len) if candidates else None

    def _from_zip() -> bytes | None:
        with zipfile.ZipFile(archive) as zip_file:
            member = _shallowest(zip_file.namelist())
            if member is None:
                return None
            with zip_file.open(member) as handle:
                return handle.read(_MAX_README_BYTES)

    def _from_tar() -> bytes | None:
        with tarfile.open(archive, "r:*") as tar_file:
            member = _shallowest(tar_file.getnames())
            if member is None:
                return None
            handle = tar_file.extractfile(member)
            if handle is None:
                return None
            with handle:
                return handle.read(_MAX_README_BYTES)

    try:
        raw = _from_zip() if archive.name.endswith(".zip") else _from_tar()
    except (OSError, tarfile.TarError, zipfile.BadZipFile, EOFError) as exc:
        logger.warning("Could not read README from %s (%s)", archive, exc)
        return None

    return None if raw is None else raw.decode("utf-8", errors="replace")


def read_original_search(archive: Path) -> dict | None:
    """The POST body that produced this archive, read back out of its README."""
    text = _readme_text(archive)
    if text is None:
        return None

    match = _README_JSON.search(text)
    if match is None:
        logger.warning("No download command JSON in the README of %s", archive)
        return None

    try:
        original_search = json.loads(match.group(1))
    except json.JSONDecodeError as exc:
        logger.warning("Unparsable download command JSON in %s (%s)", archive, exc)
        return None

    if not isinstance(original_search, dict):
        return None
    return original_search


def _download_params(original_search: dict):
    """``(protein_params, other_params, static_link)`` for a recovered search.

    Routed through the view's own serializer so a search that predates a flag
    picks up the same default the view would give it today.
    """
    serializer = DownloadStructuresSerializer(data=original_search)
    if not serializer.is_valid():
        logger.warning("Recovered search did not validate: %s", serializer.errors)
        return None, None, False
    return get_download_params(serializer.validated_data)


def _resolve_target(title: str | None, tas: str | None) -> Target | None:
    """The Target this download was built from, if it is still unambiguous.

    Titles are not unique across projects, so a title alone is only accepted
    when it matches exactly one Target.
    """
    if not title:
        return None
    targets = Target.objects.filter(title=title)
    if tas:
        targets = targets.filter(project__title=tas)
    found = list(targets[:2])
    return found[0] if len(found) == 1 else None


def _title_from_filename(file_url: str) -> str:
    """The target title an archive is named after ("A71EV2A.tar.gz" -> "A71EV2A")."""
    for suffix in _ARCHIVE_SUFFIXES:
        if file_url.endswith(suffix):
            return file_url[: -len(suffix)]
    return file_url


def _describe(directory: Path, archive: Path) -> Orphan:
    stat = archive.stat()
    create_date = datetime.fromtimestamp(stat.st_mtime, timezone.utc)
    original_search = read_original_search(archive)

    if original_search is None:
        # No readable search: the archive still names its target, which is
        # enough to restore the link itself (with --unparsed-as-static).
        return Orphan(
            directory=directory,
            archive=archive,
            task_id=directory.name,
            file_url=archive.name,
            size_bytes=stat.st_size,
            create_date=create_date,
            original_search=None,
            static_link=False,
            target=_resolve_target(_title_from_filename(archive.name), None),
            note="no download command JSON in the archive",
        )

    _, _, static_link = _download_params(original_search)
    return Orphan(
        directory=directory,
        archive=archive,
        task_id=directory.name,
        file_url=archive.name,
        size_bytes=stat.st_size,
        create_date=create_date,
        original_search=original_search,
        static_link=bool(static_link),
        target=_resolve_target(
            original_search.get("target_name"),
            original_search.get("target_access_string"),
        ),
    )


def scan_downloads() -> RestoreReport:
    """Find the download directories that no ``DownloadLinks`` row points at."""
    root = downloads_root()
    report = RestoreReport(root=root)

    if not root.is_dir():
        report.warnings.append(f"{root} does not exist")
        return report

    known = set(
        DownloadLinks.objects.exclude(task_id__isnull=True).values_list(
            "task_id", flat=True
        )
    )

    for entry in sorted(root.iterdir()):
        report.scanned += 1
        if entry.is_symlink() or not entry.is_dir():
            report.skipped.append((entry, "not a directory"))
            continue
        if entry.name in known:
            report.referenced += 1
            continue

        archive, reason = _find_archive(entry)
        if archive is None:
            report.skipped.append((entry, reason or "no archive"))
            continue

        orphan = _describe(entry, archive)
        if orphan.target is None:
            report.warnings.append(
                f"{orphan.task_id}: target could not be resolved"
                " - restoring with no target"
            )
        report.orphans.append(orphan)

    return report


def restore(
    orphans: list[Orphan],
    *,
    include_dynamic: bool = False,
    unparsed_as_static: bool = False,
) -> RestoreResult:
    """Re-create a ``DownloadLinks`` row for each orphan.

    Static links are restored with no ``keep_zip_until``, exactly as
    ``create_download()`` leaves them once they are marked static: the three
    housekeeping jobs all filter ``static_link=False``, so the archive stays.

    ``include_dynamic`` also restores dynamic links, with their original
    retention (``create_date + KEEP_UNTIL_DURATION``). That is long past, so the
    housekeeping will expire them and, after the hard-expiry grace period,
    **delete the archives**.

    ``unparsed_as_static`` restores archives with no readable search as static
    links, so the URL works again even though the search behind it is lost.
    """
    result = RestoreResult()

    for orphan in orphans:
        if not orphan.parsed and not unparsed_as_static:
            continue
        if orphan.parsed and not orphan.static_link and not include_dynamic:
            continue

        static_link = True if not orphan.parsed else orphan.static_link
        protein_params, other_params = None, None
        if orphan.original_search is not None:
            protein_params, other_params, _ = _download_params(orphan.original_search)

        try:
            with transaction.atomic():
                # Re-check inside the transaction: a concurrent download could
                # have claimed this task_id since the scan.
                if DownloadLinks.objects.filter(task_id=orphan.task_id).exists():
                    result.failed.append((orphan, "a record now exists for this task"))
                    continue
                DownloadLinks.objects.create(
                    task_id=orphan.task_id,
                    file_url=orphan.file_url,
                    static_link=static_link,
                    target=orphan.target,
                    user=None,
                    proteins=None,
                    protein_params=protein_params,
                    other_params=other_params,
                    original_search=orphan.original_search,
                    create_date=orphan.create_date,
                    keep_zip_until=(
                        None
                        if static_link
                        else orphan.create_date + KEEP_UNTIL_DURATION
                    ),
                )
        except Exception as exc:  # pylint: disable=broad-except
            # One unrestorable archive must not abandon the rest of the run.
            logger.warning("Could not restore %s (%s)", orphan.directory, exc)
            result.failed.append((orphan, str(exc)))
            continue

        logger.info(
            "Restored DownloadLinks task=%s file_url=%s static=%s target=%s",
            orphan.task_id,
            orphan.file_url,
            static_link,
            orphan.target,
        )
        result.restored.append(orphan)

    return result


def directory_size(path: Path) -> int:
    """Total bytes under ``path`` (used for reporting only)."""
    total = 0
    for dirpath, _, filenames in os.walk(path):
        for name in filenames:
            try:
                total += os.path.getsize(os.path.join(dirpath, name))
            except OSError:
                continue
    return total
