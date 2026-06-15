"""Fully delete a Target from the database and the media directory.

This is the reusable core of the target-deletion feature. It is deliberately
free of any HTTP/request concerns so it can be called from a view, a future
management command, or anywhere else in the code. The *caller* is responsible
for authorization (project membership) and for the deployment-mode guard - this
function simply does the deletion once those checks have passed.
"""

import logging
import shutil
from pathlib import Path

from django.conf import settings
from django.db import transaction
from django.db.models.deletion import ProtectedError

from viewer.models import ComputedMolecule, Target

logger = logging.getLogger(__name__)


def _computed_set_media_files(target: Target) -> list[Path]:
    """Collect the absolute paths of every computed-set file belonging to this
    target. These live under ``media/computed_set_data/`` which is shared across
    targets, so they must be removed individually - never by nuking the dir.

    Must be called *before* the database rows are deleted.
    """
    paths: list[Path] = []
    media_root = Path(settings.MEDIA_ROOT)

    for computed_set in target.computedset_set.all():
        # The written ComputedSet SDF - an absolute path stored in a TextField.
        if computed_set.written_sdf_filename:
            paths.append(Path(computed_set.written_sdf_filename))
        # The originally submitted SDF - a FileField under computed_set_data/.
        if computed_set.submitted_sdf:
            paths.append(media_root.joinpath(computed_set.submitted_sdf.name))
        # Each computed molecule's pdb file - a FileField under computed_set_data/.
        for comp_mol in computed_set.computed_molecules.all():
            if comp_mol.pdb_info:
                paths.append(media_root.joinpath(comp_mol.pdb_info.name))

    return paths


def _target_loader_media_dir(target: Target) -> Path | None:
    """The target-loader subdirectory for this target (deleted in its entirety),
    or ``None`` if the target has no ``zip_archive`` set.

    ``zip_archive`` stores the subdirectory name (e.g. ``Mpro_lb32627-66``) under
    ``media/target_loader_data/`` - see ``viewer/target_loader.py``.
    """
    subdir = str(target.zip_archive) if target.zip_archive else ""
    # Guard against a blank value so we never rmtree the parent directory.
    if not subdir.strip():
        return None
    return (
        Path(settings.MEDIA_ROOT)
        .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
        .joinpath(subdir)
    )


def delete_target(target: Target) -> None:
    """Fully delete a Target: its database graph and its media files.

    Removes the data loaded through the target loader and through the computed-set
    upload, for this target alone - all other targets are left untouched. The
    caller is responsible for authorization and the deployment-mode check.

    :param target: the Target to delete.
    """
    target_pk = target.pk
    target_title = target.title
    logger.info("Deleting target pk=%s title=%s", target_pk, target_title)

    # Collect media paths up front, while the DB rows still exist. Files are only
    # removed *after* the database transaction commits, so a DB failure can't
    # leave us with deleted files but live rows.
    computed_set_files = _computed_set_media_files(target)
    target_loader_dir = _target_loader_media_dir(target)

    with transaction.atomic():
        # 1. Computed molecules whose pdb (SiteObservation) belongs to this
        #    target. ComputedMolecule.pdb is on_delete=PROTECT, so deleting the
        #    target would raise ProtectedError if these still referenced its site
        #    observations. Remove them first.
        ComputedMolecule.objects.filter(
            pdb__cmpd__experiment__experiment_upload__target=target
        ).delete()

        # 2. Computed molecules belonging to this target's computed sets. These
        #    are linked by an m2m and so are not cleaned up by the cascade from
        #    target.delete().
        ComputedMolecule.objects.filter(computed_set__target=target).delete()

        # 3. The target itself. This cascades to the experiment uploads,
        #    experiments, site observations, computed sets, session projects,
        #    tags, result uploads, plot data, download links and jobs.
        try:
            target.delete()
        except ProtectedError:
            # A relation we don't yet handle is still protecting the target.
            # Surface which objects blocked it so the missing relation can be
            # added here, rather than failing opaquely.
            logger.exception(
                "ProtectedError deleting target pk=%s title=%s", target_pk, target_title
            )
            raise

    # DB rows are gone - now remove the media files.
    if target_loader_dir is not None:
        logger.info("Removing target loader directory %s", target_loader_dir)
        shutil.rmtree(target_loader_dir, ignore_errors=True)

    for file_path in computed_set_files:
        try:
            if file_path.is_file():
                file_path.unlink()
        except OSError:
            # A missing or unremovable computed-set file shouldn't fail the whole
            # deletion - the DB rows are already gone.
            logger.warning("Could not remove computed set file %s", file_path)

    logger.info("Deleted target pk=%s title=%s", target_pk, target_title)
