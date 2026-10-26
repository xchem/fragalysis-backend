# pylint: disable=redefined-outer-name
"""Tests for the S3Storage Django storage class."""
import datetime

import pytest
from botocore.exceptions import ClientError
from django.core.exceptions import ImproperlyConfigured, SuspiciousFileOperation
from django.core.files.base import ContentFile

from frag_s3_storage import S3Storage

from .conftest import BUCKET, REGION


def test_missing_bucket_name_is_an_error():
    with pytest.raises(ImproperlyConfigured):
        S3Storage()


def test_save_then_open_returns_the_same_content(storage):
    name = storage.save("dir/a.bin", ContentFile(b"\x00\x01binary"))

    assert name == "dir/a.bin"
    with storage.open(name) as f:
        assert f.read() == b"\x00\x01binary"


def test_save_text_content_is_stored_as_utf8(storage):
    storage.save("a.txt", ContentFile("café"))

    with storage.open("a.txt") as f:
        assert f.read() == "café".encode("utf-8")


def test_save_sets_the_content_type_from_the_name(storage, s3_client):
    storage.save("a.txt", ContentFile(b"x"))

    head = s3_client.head_object(Bucket=BUCKET, Key="a.txt")
    assert head["ContentType"] == "text/plain"


def test_save_does_not_overwrite_an_existing_file(storage):
    first = storage.save("a.txt", ContentFile(b"one"))
    second = storage.save("a.txt", ContentFile(b"two"))

    assert first == "a.txt"
    assert second != first
    with storage.open(first) as f:
        assert f.read() == b"one"


def test_open_a_missing_file_raises_file_not_found(storage):
    with pytest.raises(FileNotFoundError):
        storage.open("missing.txt")


def test_open_for_writing_is_not_supported(storage):
    with pytest.raises(ValueError):
        storage.open("a.txt", "wb")


def test_exists(storage):
    assert not storage.exists("a.txt")
    storage.save("a.txt", ContentFile(b"x"))
    assert storage.exists("a.txt")


def test_delete_removes_the_file(storage):
    storage.save("a.txt", ContentFile(b"x"))

    storage.delete("a.txt")

    assert not storage.exists("a.txt")


def test_delete_a_missing_file_is_not_an_error(storage):
    storage.delete("missing.txt")


def test_size(storage):
    storage.save("a.txt", ContentFile(b"12345"))

    assert storage.size("a.txt") == 5


def test_size_of_a_missing_file_raises_file_not_found(storage):
    with pytest.raises(FileNotFoundError):
        storage.size("missing.txt")


def test_get_modified_time_is_timezone_aware(storage):
    storage.save("a.txt", ContentFile(b"x"))

    modified = storage.get_modified_time("a.txt")

    assert isinstance(modified, datetime.datetime)
    assert modified.tzinfo is not None


def test_listdir_returns_directories_and_files(storage):
    storage.save("top.txt", ContentFile(b"x"))
    storage.save("d1/a.txt", ContentFile(b"x"))
    storage.save("d1/sub/b.txt", ContentFile(b"x"))
    storage.save("d2/c.txt", ContentFile(b"x"))

    assert storage.listdir("") == (["d1", "d2"], ["top.txt"])
    assert storage.listdir("d1") == (["sub"], ["a.txt"])
    assert storage.listdir("d1/") == (["sub"], ["a.txt"])


def test_url_is_a_presigned_get_url(storage):
    storage.save("dir/a.txt", ContentFile(b"x"))

    url = storage.url("dir/a.txt")

    assert BUCKET in url
    assert "dir/a.txt" in url
    assert "X-Amz-Signature=" in url or "Signature=" in url


def test_url_expiry_is_configurable(s3_client):  # pylint: disable=unused-argument
    storage = S3Storage(bucket_name=BUCKET, region_name=REGION, url_expiry=60)

    assert "X-Amz-Expires=60" in storage.url("a.txt")


def test_path_is_not_supported(storage):
    with pytest.raises(NotImplementedError):
        storage.path("a.txt")


def test_location_prefixes_every_key(s3_client):
    storage = S3Storage(bucket_name=BUCKET, region_name=REGION, location="media")

    storage.save("a.txt", ContentFile(b"x"))

    keys = [o["Key"] for o in s3_client.list_objects_v2(Bucket=BUCKET)["Contents"]]
    assert keys == ["media/a.txt"]
    assert storage.exists("a.txt")
    assert storage.listdir("") == ([], ["a.txt"])


@pytest.mark.parametrize("name", ["../a.txt", "dir/../../a.txt", "/abs/a.txt"])
def test_names_outside_the_location_are_rejected(storage, name):
    with pytest.raises(SuspiciousFileOperation):
        storage.exists(name)


def test_deconstruct_supports_migrations():
    storage = S3Storage(bucket_name=BUCKET, location="media")

    path, args, kwargs = storage.deconstruct()

    assert path == "frag_s3_storage.S3Storage"
    assert not args
    assert kwargs == {"bucket_name": BUCKET, "location": "media"}


def test_configured_through_the_storages_setting(settings, s3_client):
    # pylint: disable=unused-argument
    from django.core.files.storage import storages

    settings.STORAGES = {
        **settings.STORAGES,
        "default": {
            "BACKEND": "frag_s3_storage.S3Storage",
            "OPTIONS": {"bucket_name": BUCKET, "region_name": REGION},
        },
    }

    assert isinstance(storages["default"], S3Storage)
    assert storages["default"].bucket_name == BUCKET


def _access_denied(*args, **kwargs):
    raise ClientError(
        {"Error": {"Code": "AccessDenied", "Message": "Access Denied"}}, "operation"
    )


@pytest.mark.parametrize(
    "method, client_call",
    [
        ("exists", "head_object"),
        ("size", "head_object"),
        ("open", "download_fileobj"),
    ],
)
def test_errors_other_than_not_found_are_raised(storage, mocker, method, client_call):
    mocker.patch.object(storage.client, client_call, side_effect=_access_denied)

    with pytest.raises(ClientError):
        getattr(storage, method)("a.txt")
