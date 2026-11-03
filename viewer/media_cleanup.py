"""Find (and remove) orphaned debris in the media storage (see viewer.media_storage).

This is the reusable core of the ``cleanup_media`` management command, kept free
of any CLI concerns so the same reference-resolution logic can be shared with
``viewer.target_delete`` - the deleter and the reclaimer must never disagree
about what counts as a live reference.

Two areas are covered, each with an exact database cross-reference:

``target_loader_data/``
    Every top-level *directory* belongs to a target and its name is stored in
    ``Target.zip_archive`` (a FileField misused to hold a bare directory name,
    e.g. ``A71EV2A_lb32627-66`` - see ``viewer.target_loader``). A directory no
    ``Target`` names is debris, left behind by the targets that used to be
    deleted by hand in the Django shell.

``computed_set_data/``
    A flat, target-shared directory of computed-set files. A file no row points
    at is debris.

Resolving a computed-set reference is the subtle part: the same column holds
paths in three different flavours.

* absolute - a legacy ``ComputedSet.written_sdf_filename`` (a TextField; migration
  0174 made the values storage names);
* ``computed_set_data/<file>`` - ``SiteObservation.virtual_pdb_info``, and the
  ``ComputedSet.submitted_sdf`` of an assay-data upload;
* a bare basename - the ``ComputedSet.submitted_sdf`` written by
  ``viewer.cset_upload``, and the ``SiteObservation.virtual_ligand_mol`` of
  every row created by migration ``0149``, whose files *are* in
  ``computed_set_data/``.

That last flavour is why the referencing fields are listed explicitly below
rather than discovered by walking model metadata: resolving a bare name as a
storage name lands outside ``computed_set_data/`` and would mark every migrated
``.mol`` file as debris. ``unexpected_computed_set_fields`` is the guard against
the list going stale - a test asserts it stays empty.

Each area is listed once, and files are removed in batches (S3 requests are slow).
``expire_tmp_files`` removes old objects in ``tmp/``.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import PurePosixPath

from django.apps import apps
from django.conf import settings
from django.core.exceptions import SuspiciousFileOperation
from django.db import models
from django.db.models import Q

from viewer import media_storage
from viewer.models import ComputedSet, Experiment, SiteObservation, Target

logger = logging.getLogger(__name__)

# Every field that can point at a file in computed_set_data/. Keep in step with
# viewer/models.py - unexpected_computed_set_fields() fails the test-suite if a
# new FileField targets that directory without being added here.
_COMPUTED_SET_FIELDS = (
    (ComputedSet, "submitted_sdf"),
    (ComputedSet, "written_sdf_filename"),
    (SiteObservation, "virtual_pdb_info"),
    (SiteObservation, "virtual_ligand_mol"),
)

# A sample of the fields that point *inside* a target-loader directory. Used
# only as a sanity check on a directory that already looks orphaned, so a hit
# means the database is inconsistent - not that the directory is in use.
_TARGET_LOADER_PROBE_FIELDS = (
    "bound_file",
    "apo_file",
    "ligand_mol",
    "ligand_sdf",
    "virtual_ligand_mol",
)

# Legacy rows carry these in a file field in place of a real path. See
# viewer/management/commands/restore_v3_files_from_v2_archive.py.
_EMPTY_VALUES = frozenset(("", "None"))


@dataclass(frozen=True)
class Candidate:
    """One orphaned entry, ready to be removed."""

    name: str  # The storage name (i.e. "target_loader_data/<dir>")
    is_dir: bool
    size_bytes: int
    # Newest modification time (a timestamp) of the entry's files
    mtime: float


@dataclass
class AreaReport:
    """What a scan of one media subdirectory found."""

    area: str
    root: str  # The storage name of the area
    scanned: int = 0
    referenced: int = 0
    candidates: list[Candidate] = field(default_factory=list)
    # (name, reason) - entries deliberately left alone, e.g. the file watcher's
    # log, recently changed entries, anything of an unexpected type.
    skipped: list[tuple[str, str]] = field(default_factory=list)
    # Referenced names with no file in the storage. The inverse orphan; free to
    # collect here and worth knowing about.
    missing: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def reclaimable_bytes(self) -> int:
        return sum(c.size_bytes for c in self.candidates)


@dataclass
class _Entry:
    """A top-level entry of an area: a file, or a "directory" of files."""

    is_dir: bool = False
    size_bytes: int = 0
    mtime: float = 0.0


def _area_entries(area: str) -> dict[str, _Entry]:
    """The area's top-level entries by name, from one listing of the storage."""
    entries: dict[str, _Entry] = {}
    for stored in media_storage.list_prefix(area):
        first, _, rest = stored.name[len(area) + 1 :].partition("/")
        entry = entries.setdefault(first, _Entry())
        entry.is_dir = entry.is_dir or bool(rest)
        entry.size_bytes += stored.size
        entry.mtime = max(entry.mtime, stored.modified.timestamp())
    return dict(sorted(entries.items()))


def _computed_set_name(value) -> str | None:
    """The name of ``value`` within ``computed_set_data/``, or None.

    Accepts all three storage flavours described in the module docstring, and
    returns the *first* path component so a value naming something in a
    subdirectory protects that whole subdirectory (the directory is flat today,
    but this way it stays correct if it ever isn't).
    """
    name = str(value or "").strip()
    if name in _EMPTY_VALUES:
        return None

    if PurePosixPath(name).is_absolute():
        # A legacy absolute path: its storage name, if it is under MEDIA_ROOT
        try:
            name = media_storage.to_name(name)
        except SuspiciousFileOperation:
            return None  # Not a media file

    pure = PurePosixPath(name)
    if pure.parent == PurePosixPath("."):
        # A bare basename means computed_set_data/<name>.
        return pure.name

    # A storage name: either computed_set_data/... (ours) or
    # target_loader_data/... (a current-upload virtual_ligand_mol, which
    # protects nothing here).
    try:
        relative = pure.relative_to(settings.COMPUTED_SET_MEDIA_DIRECTORY)
    except ValueError:
        return None
    return relative.parts[0] if relative.parts else None


def computed_set_file_names(computed_sets, site_observations) -> set[str]:
    """Names within ``computed_set_data/`` referenced by the given querysets.

    Shared with ``viewer.target_delete``, which passes one target's rows where
    the cleanup scan passes every row.
    """
    names: set[str] = set()

    for name, ordinal, submitted, written in computed_sets.values_list(
        "name", "md_ordinal", "submitted_sdf", "written_sdf_filename"
    ):
        for value in (submitted, written):
            if resolved := _computed_set_name(value):
                names.add(resolved)

        # cset_upload renames the uploaded SDF to this and stores only the
        # *original* basename in submitted_sdf, so reconstruct the real name
        # rather than relying on written_sdf_filename being populated. Only the
        # current basename is protected: an earlier upload under the same
        # (name, ordinal) genuinely is debris.
        if submitted and name and ordinal is not None:
            names.add(f"{name}_upload_{ordinal}_{PurePosixPath(str(submitted)).name}")

    for pdb_info, ligand_mol in site_observations.values_list(
        "virtual_pdb_info", "virtual_ligand_mol"
    ):
        for value in (pdb_info, ligand_mol):
            if resolved := _computed_set_name(value):
                names.add(resolved)

    return names


def referenced_computed_set_names() -> set[str]:
    """Every name in ``computed_set_data/`` the database still points at."""
    return computed_set_file_names(
        ComputedSet.objects.all(), SiteObservation.objects.all()
    )


def referenced_target_dir_names() -> tuple[set[str], list[str]]:
    """Directory names in ``target_loader_data/`` that a Target still claims.

    Returns ``(names, warnings)``. A zip_archive that isn't a plain directory
    name is never used for matching - it goes into the warnings instead, so a
    hand-edited value can't be turned into a path traversal.
    """
    names: set[str] = set()
    warnings: list[str] = []

    for pk, value in (
        Target.objects.exclude(zip_archive__isnull=True)
        .exclude(zip_archive="")
        .values_list("pk", "zip_archive")
    ):
        candidate = str(value or "").strip()
        if candidate in _EMPTY_VALUES:
            warnings.append(f"Target pk={pk}: empty zip_archive ({value!r}) - ignored")
        elif "/" in candidate or candidate.startswith("."):
            warnings.append(f"Target pk={pk}: unusable zip_archive {value!r} - ignored")
        else:
            names.add(candidate)

    return names, warnings


def unexpected_computed_set_fields() -> list[str]:
    """FileFields pointing at ``computed_set_data/`` that aren't accounted for.

    The tripwire for ``_COMPUTED_SET_FIELDS`` going stale. Historical models are
    skipped - they mirror their originals' fields.
    """
    known = {(model.__name__, name) for model, name in _COMPUTED_SET_FIELDS}
    unexpected = []

    for model in apps.get_models():
        meta = model._meta  # pylint: disable=protected-access
        if meta.model_name.startswith("historical"):
            continue
        for model_field in meta.get_fields():
            if not isinstance(model_field, models.FileField):
                continue
            upload_to = str(getattr(model_field, "upload_to", ""))
            if settings.COMPUTED_SET_MEDIA_DIRECTORY not in upload_to:
                continue
            if (model.__name__, model_field.name) not in known:
                unexpected.append(f"{model.__name__}.{model_field.name}")

    return sorted(unexpected)


def _rows_reference_directory(dir_name: str) -> bool:
    """True if any row names a file inside this target-loader directory.

    Only ever called for a directory that already looks orphaned, so a hit is a
    database inconsistency worth reporting rather than deleting.
    """
    prefix = f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/{dir_name}/"
    query = Q()
    for name in _TARGET_LOADER_PROBE_FIELDS:
        query |= Q(**{f"{name}__startswith": prefix})
    return (
        SiteObservation.objects.filter(query).exists()
        or Experiment.objects.filter(pdb_info__startswith=prefix).exists()
    )


def _classify(name: str, entry: _Entry, want_dir: bool, newest: float) -> str | None:
    """Why this entry can't be a candidate, or None if it can be.

    An entry of the wrong type is reported rather than touched. In
    target_loader_data that is what protects the file watcher's
    ``deletions-<host>.log``. A recently changed entry can belong to an upload
    that is still in progress (see scan_target_loader_data()).
    """
    if name.startswith("."):
        return "hidden"
    if entry.is_dir != want_dir:
        return "not a directory" if want_dir else "unexpected directory"
    if entry.mtime > newest:
        return "recently changed"
    return None


def _newest_allowed(min_age: timedelta) -> float:
    return (datetime.now(timezone.utc) - min_age).timestamp()


def scan_target_loader_data(min_age: timedelta = timedelta(0)) -> AreaReport:
    """Find target-loader directories no Target claims.

    The loader stores an upload's files as the last step of its transaction, so
    for a new target the directory exists shortly before the Target that claims
    it is committed. An entry changed within ``min_age`` is therefore skipped:
    the cleanup_media command uses a safe default.
    """
    area = settings.TARGET_LOADER_MEDIA_DIRECTORY
    referenced, warnings = referenced_target_dir_names()
    report = AreaReport(area=area, root=area, warnings=warnings)
    newest = _newest_allowed(min_age)

    found = set()
    for entry_name, entry in _area_entries(area).items():
        report.scanned += 1
        name = f"{area}/{entry_name}"

        if entry_name in referenced:
            found.add(entry_name)
            report.referenced += 1
            continue

        if reason := _classify(entry_name, entry, want_dir=True, newest=newest):
            report.skipped.append((name, reason))
            continue

        if _rows_reference_directory(entry_name):
            report.skipped.append((name, "rows reference files inside it"))
            report.warnings.append(
                f"{entry_name}: no Target claims this directory but rows still"
                " reference files inside it - left alone"
            )
            continue

        report.candidates.append(
            Candidate(
                name=name, is_dir=True, size_bytes=entry.size_bytes, mtime=entry.mtime
            )
        )

    report.missing = sorted(referenced - found)
    return report


def scan_computed_set_data(min_age: timedelta = timedelta(0)) -> AreaReport:
    """Find computed-set files no row points at.

    A file changed within ``min_age`` is skipped (it can belong to a computed-set
    upload that is still in progress).
    """
    area = settings.COMPUTED_SET_MEDIA_DIRECTORY
    referenced = referenced_computed_set_names()
    report = AreaReport(area=area, root=area)
    newest = _newest_allowed(min_age)

    found = set()
    for entry_name, entry in _area_entries(area).items():
        report.scanned += 1
        name = f"{area}/{entry_name}"

        if entry_name in referenced:
            found.add(entry_name)
            report.referenced += 1
            continue

        if reason := _classify(entry_name, entry, want_dir=False, newest=newest):
            report.skipped.append((name, reason))
            continue

        report.candidates.append(
            Candidate(
                name=name, is_dir=False, size_bytes=entry.size_bytes, mtime=entry.mtime
            )
        )

    report.missing = sorted(referenced - found)
    return report


def _referenced_for(area: str) -> set[str]:
    """The live reference set for an area, re-queried."""
    if area == settings.TARGET_LOADER_MEDIA_DIRECTORY:
        return referenced_target_dir_names()[0]
    return referenced_computed_set_names()


def delete_candidates(report: AreaReport) -> tuple[int, int, list[str]]:
    """Remove a report's candidates. Returns ``(deleted, failed, errors)``.

    Every safety check is repeated here rather than trusted from scan time: an
    upload can commit between the scan and the deletion, so the reference set is
    re-queried and each entry is listed again before anything is removed.
    Failures are collected, never raised - a half-removed directory simply turns
    up as a candidate again on the next run. Files are removed in batches.
    """
    referenced = _referenced_for(report.area)
    entries = _area_entries(report.area)
    deleted = 0
    failed = 0
    errors: list[str] = []
    files: list[str] = []

    for candidate in report.candidates:
        entry_name = PurePosixPath(candidate.name).name
        entry = entries.get(entry_name)
        if (
            entry is None
            or entry.is_dir != candidate.is_dir
            or entry_name in referenced
        ):
            failed += 1
            errors.append(f"{candidate.name}: changed since the scan - left alone")
            continue

        # Logged before the event: filewatcher.sh exists because files have
        # gone missing unattributably, and this command must not become a new
        # suspect.
        logger.info(
            "cleanup_media removing %s (%s bytes)", candidate.name, candidate.size_bytes
        )
        if not candidate.is_dir:
            files.append(candidate.name)
            continue
        try:
            media_storage.delete_prefix(candidate.name)
        except OSError as exc:
            failed += 1
            errors.append(f"{candidate.name}: {exc}")
            continue
        deleted += 1

    if files:
        try:
            media_storage.delete_many(files)
            deleted += len(files)
        except OSError as exc:
            failed += len(files)
            errors.append(f"{len(files)} files: {exc}")

    return deleted, failed, errors


def expire_tmp_files(now: datetime | None = None) -> int:
    """Deletes the tmp/ objects older than MEDIA_TMP_MAX_AGE_H hours.

    tmp/ holds staged uploads, computed-set and Squonk job files, and what a
    failed upload left behind. Run periodically (see start_download_cleanup).
    Returns the number of objects deleted.
    """
    now = now or datetime.now(timezone.utc)
    cutoff = now - timedelta(hours=settings.MEDIA_TMP_MAX_AGE_H)
    expired = [o.name for o in media_storage.list_prefix("tmp") if o.modified < cutoff]
    media_storage.delete_many(expired)
    if expired:
        logger.info("Deleted %s tmp/ objects older than %s", len(expired), cutoff)
    return len(expired)
