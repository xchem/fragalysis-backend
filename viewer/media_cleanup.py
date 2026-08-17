"""Find (and remove) orphaned debris in the media directory.

This is the reusable core of the ``cleanup_media`` management command, kept free
of any CLI concerns so the same reference-resolution logic can be shared with
``viewer.target_delete`` - the deleter and the reclaimer must never disagree
about what counts as a live reference.

Two areas are covered, each with an exact database cross-reference:

``media/target_loader_data/``
    Every top-level *directory* belongs to a target and its name is stored in
    ``Target.zip_archive`` (a FileField misused to hold a bare directory name,
    e.g. ``A71EV2A_lb32627-66`` - see ``viewer.target_loader``). A directory no
    ``Target`` names is debris, left behind by the targets that used to be
    deleted by hand in the Django shell.

``media/computed_set_data/``
    A flat, target-shared directory of computed-set files. A file no row points
    at is debris.

Resolving a computed-set reference is the subtle part: the same column holds
paths in three different flavours.

* absolute - ``ComputedSet.written_sdf_filename`` (a TextField);
* ``computed_set_data/<file>`` - ``SiteObservation.virtual_pdb_info``, and the
  ``ComputedSet.submitted_sdf`` of an assay-data upload;
* a bare basename - the ``ComputedSet.submitted_sdf`` written by
  ``viewer.cset_upload``, and the ``SiteObservation.virtual_ligand_mol`` of
  every row created by migration ``0149``, whose files *are* in
  ``computed_set_data/``.

That last flavour is why the referencing fields are listed explicitly below
rather than discovered by walking model metadata: resolving a bare name against
MEDIA_ROOT lands outside ``computed_set_data/`` and would mark every migrated
``.mol`` file as debris. ``unexpected_computed_set_fields`` is the guard against
the list going stale - a test asserts it stays empty.
"""

import logging
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

from django.apps import apps
from django.conf import settings
from django.db import models
from django.db.models import Q

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

    path: Path
    is_dir: bool
    size_bytes: int
    # Reported so a human can eyeball the listing. Deliberately *not* used to
    # decide anything - see the note on in-flight uploads in AreaReport.
    mtime: float


@dataclass
class AreaReport:
    """What a scan of one media subdirectory found."""

    area: str
    root: Path
    scanned: int = 0
    referenced: int = 0
    candidates: list[Candidate] = field(default_factory=list)
    # (path, reason) - entries deliberately left alone, e.g. the file watcher's
    # log, symlinks, anything of an unexpected type.
    skipped: list[tuple[Path, str]] = field(default_factory=list)
    # Referenced names with no file on disk. The inverse orphan; free to collect
    # here and worth knowing about.
    missing: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def reclaimable_bytes(self) -> int:
        return sum(c.size_bytes for c in self.candidates)


def media_subdir(name: str) -> Path:
    """``MEDIA_ROOT/<name>``, read from settings at call time.

    Nothing in this module may cache MEDIA_ROOT at import time or the test-suite
    cannot retarget it.
    """
    return Path(settings.MEDIA_ROOT).joinpath(name)


def is_within(path: Path, root: Path) -> bool:
    """True if ``path`` resolves to somewhere inside ``root``."""
    try:
        return path.resolve().is_relative_to(root.resolve())
    except (OSError, RuntimeError):
        # Unreadable path or a symlink loop: treat as outside, i.e. hands off.
        return False


def entry_stats(path: Path, is_dir: bool) -> tuple[int, float]:
    """Total size in bytes and newest mtime of ``path``, never following symlinks.

    Symlinked files contribute their own (tiny) size, never the size of whatever
    they point at, so a staging tree of download symlinks can't inflate the
    reclaimable total.
    """
    try:
        stat = path.lstat()
    except OSError:
        return 0, 0.0

    if not is_dir:
        return stat.st_size, stat.st_mtime

    total = 0
    newest = stat.st_mtime
    for dirpath, _, filenames in os.walk(path, followlinks=False):
        for filename in filenames:
            file_path = Path(dirpath, filename)
            try:
                file_stat = file_path.lstat()
            except OSError:
                continue
            total += file_stat.st_size
            newest = max(newest, file_stat.st_mtime)
    return total, newest


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

    pure = PurePosixPath(name)
    computed_set_dir = settings.COMPUTED_SET_MEDIA_DIRECTORY

    if pure.is_absolute():
        root = PurePosixPath(str(media_subdir(computed_set_dir)))
    elif pure.parent != PurePosixPath("."):
        # Relative to MEDIA_ROOT: either computed_set_data/... (ours) or
        # target_loader_data/... (a current-upload virtual_ligand_mol, which
        # protects nothing here).
        root = PurePosixPath(computed_set_dir)
    else:
        # A bare basename means computed_set_data/<name>.
        return pure.name

    try:
        relative = pure.relative_to(root)
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


def _classify(entry, want_dir: bool) -> str | None:
    """Why this directory entry can't be a candidate, or None if it can be.

    Symlinks are never followed or removed - the download staging trees are
    built out of them - and an entry of the wrong type is reported rather than
    touched. In target_loader_data that is what protects the file watcher's
    ``deletions-<host>.log``.
    """
    if entry.name.startswith("."):
        return "hidden"
    if entry.is_symlink():
        return "symlink"
    if entry.is_dir(follow_symlinks=False) != want_dir:
        return "not a directory" if want_dir else "unexpected directory"
    return None


def _scandir(root: Path):
    """``os.scandir`` entries of ``root``, or an empty list if it isn't there."""
    try:
        with os.scandir(root) as entries:
            return sorted(entries, key=lambda e: e.name)
    except FileNotFoundError:
        return []


def scan_target_loader_data() -> AreaReport:
    """Find target-loader directories no Target claims.

    An in-flight upload is never a candidate: the loader does all its work in a
    temporary directory and only creates ``target_loader_data/<dir>/`` in
    ``_move_and_save_target_experiment``, after the transaction that saved
    ``Target.zip_archive`` has committed. The database reference therefore
    always exists before the directory does.
    """
    root = media_subdir(settings.TARGET_LOADER_MEDIA_DIRECTORY)
    referenced, warnings = referenced_target_dir_names()
    report = AreaReport(
        area=settings.TARGET_LOADER_MEDIA_DIRECTORY, root=root, warnings=warnings
    )

    found = set()
    for entry in _scandir(root):
        report.scanned += 1
        path = Path(entry.path)

        if reason := _classify(entry, want_dir=True):
            report.skipped.append((path, reason))
            continue

        found.add(entry.name)
        if entry.name in referenced:
            report.referenced += 1
            continue

        if _rows_reference_directory(entry.name):
            report.skipped.append((path, "rows reference files inside it"))
            report.warnings.append(
                f"{entry.name}: no Target claims this directory but rows still"
                " reference files inside it - left alone"
            )
            continue

        size, mtime = entry_stats(path, is_dir=True)
        report.candidates.append(
            Candidate(path=path, is_dir=True, size_bytes=size, mtime=mtime)
        )

    report.missing = sorted(referenced - found)
    return report


def scan_computed_set_data() -> AreaReport:
    """Find computed-set files no row points at."""
    root = media_subdir(settings.COMPUTED_SET_MEDIA_DIRECTORY)
    referenced = referenced_computed_set_names()
    report = AreaReport(area=settings.COMPUTED_SET_MEDIA_DIRECTORY, root=root)

    found = set()
    for entry in _scandir(root):
        report.scanned += 1
        path = Path(entry.path)

        if reason := _classify(entry, want_dir=False):
            report.skipped.append((path, reason))
            continue

        found.add(entry.name)
        if entry.name in referenced:
            report.referenced += 1
            continue

        size, mtime = entry_stats(path, is_dir=False)
        report.candidates.append(
            Candidate(path=path, is_dir=False, size_bytes=size, mtime=mtime)
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
    re-queried and each path is re-checked before anything is removed. Failures
    are collected, never raised - a half-removed directory simply turns up as a
    candidate again on the next run.
    """
    referenced = _referenced_for(report.area)
    deleted = 0
    failed = 0
    errors: list[str] = []

    for candidate in report.candidates:
        path = candidate.path
        if (
            path.is_symlink()
            or not is_within(path, report.root)
            or path.name in referenced
            or path.is_dir() != candidate.is_dir
        ):
            failed += 1
            errors.append(f"{path}: changed since the scan - left alone")
            continue

        # Logged before the event: filewatcher.sh exists because files have
        # gone missing unattributably, and this command must not become a new
        # suspect.
        logger.info("cleanup_media removing %s (%s bytes)", path, candidate.size_bytes)
        try:
            if candidate.is_dir:
                shutil.rmtree(path)
            else:
                path.unlink()
        except OSError as exc:
            failed += 1
            errors.append(f"{path}: {exc}")
            continue
        deleted += 1

    return deleted, failed, errors
