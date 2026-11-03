# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the media clean-up through the media storage seam (issue 1038).

``viewer.media_cleanup`` finds media that no row references (one listing for
each area), and removes it in batches. It also removes ``tmp/`` objects that are
older than MEDIA_TMP_MAX_AGE_H hours.

Each test runs with both media storages (the ``media_backend`` fixture).
"""
from datetime import datetime, timedelta, timezone
from io import StringIO

import pytest
from django.core.management import call_command

from viewer import media_storage
from viewer.checks import check_media_storage
from viewer.media_cleanup import (
    delete_candidates,
    expire_tmp_files,
    scan_computed_set_data,
    scan_target_loader_data,
)
from viewer.models import ComputedSet


@pytest.fixture
def target(user, make_project, make_target):
    target = make_target(make_project("lb1", members=[user]), title="Mpro")
    target.zip_archive = "Mpro_lb1"
    target.save()
    return target


def _all_names() -> list[str]:
    return sorted(o.name for o in media_storage.list_prefix(""))


# target_loader_data -------------------------------------------------------------------------


@pytest.fixture
def loader_media(media_backend, target):
    media_storage.store_bytes("target_loader_data/Mpro_lb1/upload_1/a.pdb", b"a")
    media_storage.store_bytes("target_loader_data/Orphan/upload_1/b.pdb", b"bb")
    media_storage.store_bytes("target_loader_data/Orphan/upload_1/sub/c.pdb", b"ccc")
    media_storage.store_bytes("target_loader_data/deletions-stack-0.log", b"log")
    media_storage.store_bytes("target_loader_data/.hidden/x", b"x")


def test_scan_target_loader_data_finds_unclaimed_directories(loader_media):
    report = scan_target_loader_data()

    assert [(c.name, c.is_dir, c.size_bytes) for c in report.candidates] == [
        ("target_loader_data/Orphan", True, 5)
    ]
    assert report.referenced == 1
    assert sorted(reason for _, reason in report.skipped) == [
        "hidden",
        "not a directory",
    ]
    assert report.missing == []


def test_scan_target_loader_data_reports_a_missing_directory(media_backend, target):
    report = scan_target_loader_data()

    assert report.missing == ["Mpro_lb1"]
    assert report.candidates == []


def test_delete_candidates_removes_unclaimed_directories(loader_media):
    report = scan_target_loader_data()

    deleted, failed, errors = delete_candidates(report)

    assert (deleted, failed, errors) == (1, 0, [])
    assert _all_names() == [
        "target_loader_data/.hidden/x",
        "target_loader_data/Mpro_lb1/upload_1/a.pdb",
        "target_loader_data/deletions-stack-0.log",
    ]


def test_delete_candidates_keeps_a_directory_claimed_after_the_scan(
    loader_media, make_project, make_target
):
    report = scan_target_loader_data()
    late = make_target(make_project("lb2"), title="Orphan")
    late.zip_archive = "Orphan"
    late.save()

    deleted, failed, _ = delete_candidates(report)

    assert (deleted, failed) == (0, 1)
    assert media_storage.exists("target_loader_data/Orphan/upload_1/b.pdb")


def test_the_command_leaves_recent_entries_alone(loader_media):
    out = StringIO()

    call_command("cleanup_media", "--delete", stdout=out, stderr=out)

    assert media_storage.exists("target_loader_data/Orphan/upload_1/b.pdb")
    assert "recently changed" in out.getvalue()


# computed_set_data --------------------------------------------------------------------------


@pytest.fixture
def computed_set_media(media_backend, target, user):
    ComputedSet.objects.create(
        name="cset",
        md_ordinal=1,
        target=target,
        owner_user=user,
        submitted_sdf="hits.sdf",
        written_sdf_filename="computed_set_data/cset_upload_1_hits.sdf",
    )
    media_storage.store_bytes("computed_set_data/cset_upload_1_hits.sdf", b"s")
    media_storage.store_bytes("computed_set_data/orphan.sdf", b"oo")
    media_storage.store_bytes("computed_set_data/sub/deep.sdf", b"d")


def test_scan_computed_set_data_finds_unreferenced_files(computed_set_media):
    report = scan_computed_set_data()

    assert [(c.name, c.is_dir, c.size_bytes) for c in report.candidates] == [
        ("computed_set_data/orphan.sdf", False, 2)
    ]
    assert [reason for _, reason in report.skipped] == ["unexpected directory"]


def test_delete_candidates_removes_unreferenced_files(computed_set_media):
    deleted, failed, _ = delete_candidates(scan_computed_set_data())

    assert (deleted, failed) == (1, 0)
    assert not media_storage.exists("computed_set_data/orphan.sdf")
    assert media_storage.exists("computed_set_data/cset_upload_1_hits.sdf")


# tmp/ expiry ------------------------------------------------------------------------------


def test_expire_tmp_files_deletes_only_old_objects(media_backend, settings):
    settings.MEDIA_TMP_MAX_AGE_H = 48
    media_storage.store_bytes("tmp/uploads/u1/Mpro.tgz", b"old")
    media_storage.store_bytes("tmp/cset/u2/set.sdf", b"old")
    media_storage.store_bytes("downloads/t1/Mpro.tar.gz", b"not tmp")
    now = datetime.now(timezone.utc)

    assert expire_tmp_files(now=now + timedelta(hours=47)) == 0
    assert expire_tmp_files(now=now + timedelta(hours=49)) == 2

    assert _all_names() == ["downloads/t1/Mpro.tar.gz"]


def test_expire_tmp_files_with_no_tmp_directory(media_backend):
    assert expire_tmp_files() == 0


# The setting -------------------------------------------------------------------------------------


def test_a_too_short_tmp_age_fails_the_start_up_check(settings):
    settings.MEDIA_STORAGE = "FILESYSTEM"
    settings.MEDIA_TMP_MAX_AGE_H = 0

    assert [e.id for e in check_media_storage(None)] == ["viewer.E004"]


def test_the_default_tmp_age_passes_the_start_up_check(settings):
    settings.MEDIA_STORAGE = "FILESYSTEM"
    settings.MEDIA_TMP_MAX_AGE_H = 48

    assert check_media_storage(None) == []
