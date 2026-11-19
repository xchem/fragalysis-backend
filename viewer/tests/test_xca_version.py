"""The XChemAlign version recorded per upload (issue m2ms-2437).

Covers the loader's extraction of ``xca_git_info->tag``, the data migration
that back-fills already loaded uploads from their own ``meta_aligner.yaml``,
and ``Target.get_xca_version``.
"""
# pylint: disable=redefined-outer-name,unused-argument
import importlib
from datetime import datetime, timezone

import pytest
import yaml
from django.apps import apps as django_apps

from viewer.models import ExperimentUpload
from viewer.target_loader import _read_xca_version

populate_xca_version = importlib.import_module(
    "viewer.migrations.0175_populate_experimentupload_xca_version"
).populate_xca_version


@pytest.fixture
def target(make_project, make_target, user):
    project = make_project("proposal", members=[user])
    target = make_target(project, title="Xca")
    target.zip_archive = "Xca_proposal"
    target.save()
    return target


@pytest.fixture
def make_upload(target, user):
    def _make(n, xca_version=None, upload_data_dir=None):
        return ExperimentUpload.objects.create(
            project=target.project,
            target=target,
            committer=user,
            commit_datetime=datetime(2026, 1, n, tzinfo=timezone.utc),
            upload_data_dir=upload_data_dir or f"upload_{n}",
            upload_version=n,
            file=f"bundle_upload_{n}.tgz",
            xca_version=xca_version,
        )

    return _make


@pytest.fixture
def write_meta(settings, tmp_path, target):
    settings.MEDIA_ROOT = str(tmp_path / "media")
    target_dir = (
        tmp_path / "media" / settings.TARGET_LOADER_MEDIA_DIRECTORY / "Xca_proposal"
    )

    def _write(upload_dir, version_number, tag, version_dir=None):
        path = target_dir / upload_dir
        path.mkdir(parents=True)
        meta = {
            "version_number": version_number,
            "version_dir": version_dir or upload_dir,
            "xca_git_info": {"branch": "master", "tag": tag, "dirty": False},
        }
        (path / "meta_aligner.yaml").write_text(yaml.safe_dump(meta))

    return _write


@pytest.mark.parametrize(
    "meta, expected",
    [
        ({"xca_git_info": {"tag": "2.3.1", "dirty": False}}, "2.3.1"),
        # outside PRODUCTION a dirty bundle loads, and its tag is still recorded
        ({"xca_git_info": {"tag": "2.3.2", "dirty": True}}, "2.3.2"),
        ({"xca_git_info": {"tag": None, "dirty": False}}, None),
        ({"xca_git_info": {"dirty": False}}, None),
        ({}, None),
    ],
)
def test_read_xca_version(meta, expected):
    assert _read_xca_version(meta) == expected


def test_target_returns_latest_uploads_version(make_upload, target):
    make_upload(1, xca_version="2.2.12")
    make_upload(2, xca_version="2.3.1")

    assert target.get_xca_version() == "2.3.1"


def test_target_skips_uploads_without_a_version(make_upload, target):
    make_upload(1, xca_version="2.2.12")
    make_upload(2)

    assert target.get_xca_version() == "2.2.12"


def test_target_without_uploads_has_no_version(target):
    assert target.get_xca_version() is None


def test_migration_matches_each_upload_to_its_own_directory(
    make_upload, write_meta, target
):
    write_meta("upload_1", 1, "2.2.12")
    write_meta("upload_2", 2, "2.3.1")
    first = make_upload(1)
    second = make_upload(2)

    populate_xca_version(django_apps, None)

    first.refresh_from_db()
    second.refresh_from_db()
    assert first.xca_version == "2.2.12"
    assert second.xca_version == "2.3.1"
    assert target.get_xca_version() == "2.3.1"


def test_migration_skips_metadata_that_belongs_to_another_upload(
    make_upload, write_meta
):
    # The directory says upload_1 but the file inside declares itself upload_2.
    write_meta("upload_1", 2, "2.3.1", version_dir="upload_2")
    upload = make_upload(1)

    populate_xca_version(django_apps, None)

    upload.refresh_from_db()
    assert upload.xca_version is None


def test_migration_leaves_missing_files_and_untagged_bundles_null(
    make_upload, write_meta
):
    write_meta("upload_2", 2, None)
    missing = make_upload(1)
    untagged = make_upload(2)

    populate_xca_version(django_apps, None)

    missing.refresh_from_db()
    untagged.refresh_from_db()
    assert missing.xca_version is None
    assert untagged.xca_version is None


def test_migration_does_not_overwrite_a_recorded_version(make_upload, write_meta):
    write_meta("upload_1", 1, "2.2.12")
    upload = make_upload(1, xca_version="2.3.0")

    populate_xca_version(django_apps, None)

    upload.refresh_from_db()
    assert upload.xca_version == "2.3.0"
