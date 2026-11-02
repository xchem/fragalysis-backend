"""Tests for migration 0149 (merge ComputedMolecule into SiteObservation).

The migration writes a .mol file under MEDIA_ROOT for each ComputedMolecule.
A new stack has no computed molecules, and in S3 media storage it has no
MEDIA_ROOT directory (issue 1038). Then the migration must not touch the
filesystem, or "manage.py migrate" fails on a new S3 stack.

ComputedMolecule is not in the current schema, so the migration function is
called with stub models.
"""
import importlib
from unittest import mock

_MIGRATION = importlib.import_module("viewer.migrations.0149_auto_20251111_1558")


def _apps_with_no_computed_molecules():
    models: dict[str, mock.MagicMock] = {}

    def get_model(_app_label, model_name):
        return models.setdefault(model_name, mock.MagicMock(name=model_name))

    apps = mock.Mock()
    apps.get_model.side_effect = get_model
    get_model("viewer", "ComputedMolecule").objects.all.return_value = []
    get_model("viewer", "ComputedMolecule").objects.exists.return_value = False
    return apps


def test_no_media_directory_is_needed_without_computed_molecules(settings, tmp_path):
    media_root = tmp_path / "missing" / "media"
    settings.MEDIA_ROOT = str(media_root)

    _MIGRATION.Migration.merge_computedmolecule(
        _apps_with_no_computed_molecules(), None
    )

    assert not (tmp_path / "missing").exists()
