# pylint: disable=redefined-outer-name
"""Tests for the S3Storage file transfer, copy, batch delete and listing methods.

Django's Storage API transfers one file at a time through Python file objects.
These methods let a caller upload and download local files directly (with
boto3's managed, multipart transfers), copy objects inside the bucket, delete
many objects in one request, and list every object under a prefix.
"""
import datetime
import logging

import pytest
from boto3.exceptions import S3UploadFailedError
from django.core.exceptions import ImproperlyConfigured, SuspiciousFileOperation
from django.core.files.base import ContentFile

from frag_s3_storage import S3Storage

from .conftest import BUCKET, REGION


@pytest.fixture
def located(s3_client):  # pylint: disable=unused-argument
    """An S3Storage with a location, so that keys and names differ."""
    return S3Storage(bucket_name=BUCKET, region_name=REGION, location="media")


def _keys(s3_client):
    return sorted(
        o["Key"] for o in s3_client.list_objects_v2(Bucket=BUCKET).get("Contents", [])
    )


# upload_file / download_file -------------------------------------------------


def test_upload_file_stores_a_local_file_under_its_name(located, s3_client, tmp_path):
    local = tmp_path / "a.pdb"
    local.write_bytes(b"ATOM\n")

    assert located.upload_file(local, "dir/a b#1.pdb") == "dir/a b#1.pdb"

    assert _keys(s3_client) == ["media/dir/a b#1.pdb"]
    assert located.open("dir/a b#1.pdb").read() == b"ATOM\n"


def test_upload_file_overwrites_an_existing_file(located, tmp_path):
    located.save("a.txt", ContentFile(b"old"))
    local = tmp_path / "a.txt"
    local.write_bytes(b"new")

    located.upload_file(local, "a.txt")

    assert located.open("a.txt").read() == b"new"


def test_upload_file_sets_the_content_type(located, s3_client, tmp_path):
    local = tmp_path / "a.png"
    local.write_bytes(b"\x89PNG")

    located.upload_file(local, "a.png")

    head = s3_client.head_object(Bucket=BUCKET, Key="media/a.png")
    assert head["ContentType"] == "image/png"


def test_download_file_writes_the_object_to_a_local_file(located, tmp_path):
    located.save("dir/a.txt", ContentFile(b"content"))
    local = tmp_path / "out.txt"

    located.download_file("dir/a.txt", local)

    assert local.read_bytes() == b"content"


def test_download_file_of_a_missing_file_raises_file_not_found(located, tmp_path):
    with pytest.raises(FileNotFoundError):
        located.download_file("missing.txt", tmp_path / "out.txt")


# upload_file retries ---------------------------------------------------------


@pytest.fixture
def retrying(s3_client):  # pylint: disable=unused-argument
    """An S3Storage that retries an upload at once (no delay)."""
    return S3Storage(
        bucket_name=BUCKET, region_name=REGION, upload_attempts=3, retry_delay_s=0
    )


def _fail_first(monkeypatch, storage, failures):
    """Makes the client's first ``failures`` uploads fail. Returns the call log."""
    real_upload = storage.client.upload_file
    calls = []

    def flaky_upload(*args, **kwargs):
        calls.append(args)
        if len(calls) <= failures:
            raise S3UploadFailedError(
                "An error occurred (SignatureDoesNotMatch) when calling"
                " the PutObject operation: None"
            )
        return real_upload(*args, **kwargs)

    monkeypatch.setattr(storage.client, "upload_file", flaky_upload)
    return calls


def test_upload_file_retries_a_failed_upload(retrying, monkeypatch, tmp_path):
    local = tmp_path / "a.pdb"
    local.write_bytes(b"ATOM\n")
    calls = _fail_first(monkeypatch, retrying, failures=2)

    assert retrying.upload_file(local, "a.pdb") == "a.pdb"

    assert len(calls) == 3
    assert retrying.open("a.pdb").read() == b"ATOM\n"


def test_upload_file_raises_when_every_attempt_fails(retrying, monkeypatch, tmp_path):
    local = tmp_path / "a.pdb"
    local.write_bytes(b"ATOM\n")
    calls = _fail_first(monkeypatch, retrying, failures=3)

    with pytest.raises(S3UploadFailedError):
        retrying.upload_file(local, "a.pdb")

    assert len(calls) == 3
    assert not retrying.exists("a.pdb")


def test_upload_file_logs_a_warning_for_each_retry(
    retrying, monkeypatch, tmp_path, caplog
):
    local = tmp_path / "a.pdb"
    local.write_bytes(b"ATOM\n")
    _fail_first(monkeypatch, retrying, failures=2)

    with caplog.at_level(logging.WARNING, logger="frag_s3_storage.storage"):
        retrying.upload_file(local, "a.pdb")

    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 2
    assert "SignatureDoesNotMatch" in warnings[0].getMessage()
    assert "a.pdb" in warnings[0].getMessage()


def test_upload_file_does_not_retry_a_missing_local_file(
    retrying, monkeypatch, tmp_path
):
    calls = _fail_first(monkeypatch, retrying, failures=0)

    with pytest.raises(FileNotFoundError):
        retrying.upload_file(tmp_path / "missing.pdb", "a.pdb")

    assert len(calls) == 1


def test_upload_file_tries_three_times_by_default(storage):
    assert storage.upload_attempts == 3


def test_an_upload_attempts_value_below_one_is_refused(
    s3_client,
):  # pylint: disable=unused-argument
    with pytest.raises(ImproperlyConfigured):
        S3Storage(bucket_name=BUCKET, region_name=REGION, upload_attempts=0)


# copy ------------------------------------------------------------------------


def test_copy_makes_a_second_object(located, s3_client):
    located.save("a.txt", ContentFile(b"x"))

    located.copy("a.txt", "dir/b.txt")

    assert _keys(s3_client) == ["media/a.txt", "media/dir/b.txt"]
    assert located.open("dir/b.txt").read() == b"x"


def test_copy_of_a_missing_file_raises_file_not_found(located):
    with pytest.raises(FileNotFoundError):
        located.copy("missing.txt", "b.txt")


# delete_many -----------------------------------------------------------------


def test_delete_many_deletes_every_name(located, s3_client):
    for name in ("a.txt", "dir/b.txt", "dir/c.txt"):
        located.save(name, ContentFile(b"x"))

    failed = located.delete_many(["a.txt", "dir/b.txt"])

    assert failed == []
    assert _keys(s3_client) == ["media/dir/c.txt"]


def test_delete_many_uses_batches_of_1000(located, mocker):
    spy = mocker.spy(located.client, "delete_objects")

    located.delete_many([f"f{i}.txt" for i in range(2001)])

    batch_sizes = [len(c.kwargs["Delete"]["Objects"]) for c in spy.call_args_list]
    assert batch_sizes == [1000, 1000, 1]


def test_delete_many_of_no_names_makes_no_request(located, mocker):
    spy = mocker.spy(located.client, "delete_objects")

    assert located.delete_many([]) == []
    assert spy.call_count == 0


def test_delete_many_returns_the_names_that_failed(located, mocker):
    mocker.patch.object(
        located.client,
        "delete_objects",
        return_value={"Errors": [{"Key": "media/dir/b.txt", "Code": "AccessDenied"}]},
    )

    assert located.delete_many(["a.txt", "dir/b.txt"]) == ["dir/b.txt"]


def test_delete_many_rejects_names_outside_the_location(located):
    with pytest.raises(SuspiciousFileOperation):
        located.delete_many(["../a.txt"])


# iter_objects ----------------------------------------------------------------


def test_iter_objects_lists_every_object_under_a_prefix(located):
    for name in ("Mpro/a.txt", "Mpro/sub/b.txt", "Mpro_2/c.txt", "top.txt"):
        located.save(name, ContentFile(b"xyz"))

    objects = sorted(located.iter_objects("Mpro/"))

    assert [(name, size) for name, size, _ in objects] == [
        ("Mpro/a.txt", 3),
        ("Mpro/sub/b.txt", 3),
    ]
    assert all(isinstance(m, datetime.datetime) for _, _, m in objects)


def test_iter_objects_with_no_prefix_lists_everything(located):
    for name in ("a.txt", "dir/b.txt"):
        located.save(name, ContentFile(b"x"))

    assert sorted(n for n, _, _ in located.iter_objects()) == ["a.txt", "dir/b.txt"]


def test_iter_objects_pages_through_many_objects(located, s3_client):
    for i in range(1005):
        s3_client.put_object(Bucket=BUCKET, Key=f"media/many/{i:04}.txt", Body=b"")

    assert len(list(located.iter_objects("many/"))) == 1005


# max_pool_connections -------------------------------------------------------


def test_max_pool_connections_is_configurable(s3_client):
    # pylint: disable=unused-argument
    storage = S3Storage(bucket_name=BUCKET, region_name=REGION, max_pool_connections=32)

    assert storage.client.meta.config.max_pool_connections == 32


def test_max_pool_connections_defaults_to_16(storage):
    assert storage.client.meta.config.max_pool_connections == 16
