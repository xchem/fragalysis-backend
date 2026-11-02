import logging
from pathlib import Path

import yaml
from django.conf import settings
from django.db import migrations

logger = logging.getLogger(__name__)

# Kept local rather than imported from viewer.target_loader: a migration must
# not change behaviour when the application code later does.
METADATA_FILE = "meta_aligner.yaml"

try:
    _YamlLoader = yaml.CSafeLoader
except AttributeError:  # pragma: no cover - libyaml not compiled in
    _YamlLoader = yaml.SafeLoader


def populate_xca_version(apps, schema_editor):
    """Fill ``ExperimentUpload.xca_version`` from each upload's own bundle.

    Every upload unpacks into its own subdirectory,
    ``<MEDIA_ROOT>/<TARGET_LOADER_MEDIA_DIRECTORY>/<target.zip_archive>/<upload_data_dir>``,
    and the ``meta_aligner.yaml`` in it is the one that upload brought. That
    file is matched to the row by directory and, as a cross-check, by the
    ``version_dir`` and ``version_number`` the file itself declares - a row
    only gets a version from a file that agrees it belongs to that upload.

    Rows whose directory or file is missing, or whose file declares no tag,
    are left null and logged.
    """
    ExperimentUpload = apps.get_model("viewer", "ExperimentUpload")

    base = Path(settings.MEDIA_ROOT).joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
    updated = []

    for upload in (
        ExperimentUpload.objects.filter(xca_version__isnull=True)
        .select_related("target")
        .order_by("pk")
    ):
        target_dir = str(upload.target.zip_archive or "").strip()
        upload_dir = str(upload.upload_data_dir or "").strip()
        if not target_dir or not upload_dir:
            logger.warning(
                "ExperimentUpload pk=%s: no target or upload directory, skipped",
                upload.pk,
            )
            continue

        path = base.joinpath(target_dir, upload_dir, METADATA_FILE)
        if not path.is_file():
            logger.warning("ExperimentUpload pk=%s: %s missing, skipped", upload.pk, path)
            continue

        try:
            with path.open(encoding="utf-8") as f:
                meta = yaml.load(f, Loader=_YamlLoader) or {}
        except (OSError, yaml.YAMLError) as exc:
            logger.warning(
                "ExperimentUpload pk=%s: cannot read %s (%s), skipped",
                upload.pk,
                path,
                exc,
            )
            continue

        declared_dir = str(meta.get("version_dir", ""))
        try:
            declared_version = int(meta.get("version_number"))
        except (TypeError, ValueError):
            declared_version = None
        if declared_dir != upload_dir or declared_version != upload.upload_version:
            logger.warning(
                "ExperimentUpload pk=%s (%s, version %s): %s declares %s, version %s;"
                " not this upload's metadata, skipped",
                upload.pk,
                upload_dir,
                upload.upload_version,
                path,
                declared_dir,
                declared_version,
            )
            continue

        xca_git_info = meta.get("xca_git_info") or {}
        tag = xca_git_info.get("tag") if isinstance(xca_git_info, dict) else None
        if tag is None:
            logger.info("ExperimentUpload pk=%s: %s has no XCA tag", upload.pk, path)
            continue

        upload.xca_version = str(tag)
        updated.append(upload)

    if updated:
        ExperimentUpload.objects.bulk_update(updated, ["xca_version"], batch_size=500)
    logger.info("Populated xca_version on %d ExperimentUpload rows", len(updated))


class Migration(migrations.Migration):
    dependencies = [
        ("viewer", "0174_experimentupload_xca_version"),
    ]

    operations = [
        migrations.RunPython(populate_xca_version, migrations.RunPython.noop),
    ]
