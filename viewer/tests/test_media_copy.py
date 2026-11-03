# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the copy_media command (issue 1038): the data migration.

copy_media copies an existing stack's media between a filesystem directory (the
NFS media directory) and the S3 bucket of MEDIA_S3_STORAGE_OPTIONS, whatever
MEDIA_STORAGE is. Object keys are the files' relative paths, so the layout does
not change. A run copies only files that are missing or have a different size,
so it can be run again (to resume, or to copy the changes at the cutover).

The filesystem side is a temporary directory and the S3 side is a moto bucket.
"""
import os
import sys
from datetime import datetime, timedelta, timezone
from io import StringIO
from pathlib import Path

import boto3
import pytest
from django.core.management import CommandError, call_command
from moto import mock_aws

from frag_s3_storage import S3Storage
from viewer.media_copy import FS_TO_S3, Listing, plan_copy
from viewer.models import DownloadLinks, SiteObservation

_BUCKET = "frag-copy-test"
_REGION = "eu-west-2"


@pytest.fixture
def bucket(settings, monkeypatch):
    """The S3 side: a moto bucket in MEDIA_S3_STORAGE_OPTIONS."""
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    with mock_aws():
        boto3.client("s3", region_name=_REGION).create_bucket(
            Bucket=_BUCKET, CreateBucketConfiguration={"LocationConstraint": _REGION}
        )
        settings.MEDIA_S3_STORAGE_OPTIONS = {
            "bucket_name": _BUCKET,
            "region_name": _REGION,
            "location": "media",
        }
        yield S3Storage(**settings.MEDIA_S3_STORAGE_OPTIONS)


@pytest.fixture
def source(tmp_path) -> Path:
    """The filesystem side (an NFS media directory)."""
    root = tmp_path / "media"
    root.mkdir()
    return root


def _write(root: Path, name: str, data: bytes = b"x") -> None:
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def _keys(bucket: S3Storage) -> dict[str, int]:
    return {name: size for name, size, _ in bucket.iter_objects()}


def _run(*argv) -> str:
    out = StringIO()
    call_command("copy_media", *argv, stdout=out, stderr=out)
    return out.getvalue()


# Copying -----------------------------------------------------------------------------


def test_files_are_copied_with_their_names(db, bucket, source):
    _write(source, "target_loader_data/Mpro_lb1/upload_1/a b#1.pdb", b"pdb")
    _write(source, "computed_set_data/cset_upload_1_hits.sdf", b"sdf")

    _run("--source", str(source))

    assert _keys(bucket) == {
        "target_loader_data/Mpro_lb1/upload_1/a b#1.pdb": 3,
        "computed_set_data/cset_upload_1_hits.sdf": 3,
    }
    assert bucket.open("computed_set_data/cset_upload_1_hits.sdf").read() == b"sdf"


def test_temporary_and_watcher_files_are_not_copied(db, bucket, source):
    _write(source, "target_loader_data/Mpro_lb1/a.pdb")
    _write(source, "target_loader_data/deletions-stack-0.log")
    _write(source, ".filewatcher-pause.123")
    _write(source, "tmp/uploads/u1/Mpro.tgz")
    _write(source, "squonk_upload/job-1/job.sdf")

    _run("--source", str(source))

    assert list(_keys(bucket)) == ["target_loader_data/Mpro_lb1/a.pdb"]


def test_tmp_is_copied_when_asked(db, bucket, source):
    _write(source, "tmp/uploads/u1/Mpro.tgz")

    _run("--source", str(source), "--include-tmp")

    assert list(_keys(bucket)) == ["tmp/uploads/u1/Mpro.tgz"]


def test_only_downloads_that_can_be_fetched_are_copied(db, bucket, source):
    now = datetime.now(timezone.utc)
    links = {
        "static": {"static_link": True, "expired_date": now},
        "dynamic": {"static_link": False, "expired_date": None},
        "expired": {"static_link": False, "expired_date": now - timedelta(days=1)},
    }
    for task_id, fields in links.items():
        DownloadLinks.objects.create(
            task_id=task_id, file_url="Mpro.tar.gz", create_date=now, **fields
        )
        _write(source, f"downloads/{task_id}/Mpro.tar.gz")
    _write(source, "downloads/no-record/Mpro.tar.gz")
    _write(source, "downloads/dicttocsv/u1/sites.csv")

    _run("--source", str(source))

    assert sorted(_keys(bucket)) == [
        "downloads/dynamic/Mpro.tar.gz",
        "downloads/static/Mpro.tar.gz",
    ]


def test_a_second_run_copies_only_what_changed(db, bucket, source):
    _write(source, "computed_set_data/a.sdf", b"a")
    _write(source, "computed_set_data/b.sdf", b"b")
    _run("--source", str(source))
    _write(source, "computed_set_data/b.sdf", b"bigger")
    _write(source, "computed_set_data/c.sdf", b"c")

    output = _run("--source", str(source))

    assert "copied 2 files" in output
    assert _keys(bucket)["computed_set_data/b.sdf"] == len(b"bigger")


def test_dry_run_copies_nothing_and_reports_each_prefix(db, bucket, source):
    _write(source, "target_loader_data/Mpro_lb1/a.pdb", b"aaa")
    _write(source, "target_loader_data/Mpro_lb1/b.pdb", b"bb")
    _write(source, "computed_set_data/a.sdf", b"a")

    output = _run("--source", str(source), "--dry-run")

    assert _keys(bucket) == {}
    assert "target_loader_data" in output and "2 files" in output
    assert "would copy 3 files" in output


def test_symlinks_are_reported_and_skipped(db, bucket, source, tmp_path):
    outside = tmp_path / "outside.pdb"
    outside.write_text("x")
    (source / "target_loader_data").mkdir()
    (source / "target_loader_data" / "link.pdb").symlink_to(outside)

    output = _run("--source", str(source))

    assert _keys(bucket) == {}
    assert "symlink" in output


def test_names_too_long_for_an_s3_key_are_reported_and_skipped(db, bucket):
    # (A local path this long is not possible on all systems, so no file is made.)
    long_name = "target_loader_data/" + "/".join(["d" * 250] * 5) + "/x.pdb"
    source = Listing(files={long_name: 1, "computed_set_data/a.sdf": 1})

    plan = plan_copy(source, Listing(), bucket, FS_TO_S3, include_tmp=False)

    assert plan.to_copy == ["computed_set_data/a.sdf"]
    assert any("longer than 1024 bytes" in p for p in plan.problems)


@pytest.mark.skipif(sys.platform == "darwin", reason="APFS needs UTF-8 file names")
def test_names_that_are_not_utf8_are_reported_and_skipped(db, bucket, source):
    (source / "computed_set_data").mkdir()
    with open(os.fsencode(source / "computed_set_data") + b"/bad\xff.sdf", "wb") as f:
        f.write(b"x")

    output = _run("--source", str(source))

    assert _keys(bucket) == {}
    assert "not UTF-8" in output


def test_s3_to_fs_copies_back(db, bucket, source):
    bucket.save("target_loader_data/Mpro_lb1/a b#1.pdb", StringIO("pdb"))

    _run("--source", str(source), "--direction", "s3-to-fs")

    assert (source / "target_loader_data/Mpro_lb1/a b#1.pdb").read_bytes() == b"pdb"


# Verifying ------------------------------------------------------------------------------


def test_verify_passes_after_a_copy(db, bucket, source):
    _write(source, "computed_set_data/a.sdf")
    _run("--source", str(source))

    output = _run("--source", str(source), "--verify")

    assert "verified" in output.lower()


def test_verify_fails_for_a_missing_or_different_file(db, bucket, source):
    _write(source, "computed_set_data/a.sdf", b"a")
    _write(source, "computed_set_data/b.sdf", b"b")
    _run("--source", str(source))
    bucket.delete("computed_set_data/a.sdf")
    _write(source, "computed_set_data/b.sdf", b"bigger")

    with pytest.raises(CommandError, match="2 files"):
        _run("--source", str(source), "--verify")


def test_verify_copies_nothing(db, bucket, source):
    _write(source, "computed_set_data/a.sdf")

    with pytest.raises(CommandError):
        _run("--source", str(source), "--verify")

    assert _keys(bucket) == {}


# References ----------------------------------------------------------------------------------


def test_references_pass_when_every_source_file_is_copied(db, bucket, source):
    _write(source, "target_loader_data/Mpro_lb1/a.pdb")
    SiteObservation.objects.create(apo_file="target_loader_data/Mpro_lb1/a.pdb")
    # Already missing on the source: reported, not a failure
    SiteObservation.objects.create(apo_file="target_loader_data/Mpro_lb1/gone.pdb")
    _run("--source", str(source))

    output = _run("--source", str(source), "--check-references")

    assert "gone.pdb" in output
    assert "missing on the source too" in output


def test_references_fail_for_a_file_that_was_not_copied(db, bucket, source):
    _write(source, "target_loader_data/Mpro_lb1/a.pdb")
    SiteObservation.objects.create(apo_file="target_loader_data/Mpro_lb1/a.pdb")

    with pytest.raises(CommandError, match="1 referenced file"):
        _run("--source", str(source), "--check-references")
