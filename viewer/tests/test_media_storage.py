# pylint: disable=redefined-outer-name
"""Tests for the media storage configuration (issue 1038).

MEDIA_STORAGE selects where media files go: the local filesystem (the default)
or an S3 bucket. The settings build Django's STORAGES from the environment,
and a system check reports a bad configuration at start-up.
"""
import importlib

import boto3
import pytest
from django.core.checks import Error, registry
from moto import mock_aws

from viewer.checks import check_media_storage

_FILESYSTEM_BACKEND = "django.core.files.storage.FileSystemStorage"
_S3_BACKEND = "frag_s3_storage.S3Storage"
_BUCKET = "frag-media-test"
_REGION = "eu-west-2"

_MEDIA_ENV = [
    "MEDIA_STORAGE",
    "MEDIA_S3_BUCKET_NAME",
    "MEDIA_S3_LOCATION",
    "MEDIA_S3_ENDPOINT_URL",
    "MEDIA_S3_REGION_NAME",
    "MEDIA_S3_ACCESS_KEY_ID",
    "MEDIA_S3_SECRET_ACCESS_KEY",
    "MEDIA_S3_URL_EXPIRY_S",
]


def _load_settings(monkeypatch, **env):
    """Re-executes fragalysis.settings with only the given media variables set."""
    for name in _MEDIA_ENV:
        monkeypatch.delenv(name, raising=False)
    for name, value in env.items():
        monkeypatch.setenv(name, value)
    return importlib.reload(importlib.import_module("fragalysis.settings"))


@pytest.fixture
def s3_bucket(monkeypatch):
    """A moto-backed S3 bucket. Fake credentials keep real AWS out of reach."""
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    with mock_aws():
        boto3.client("s3", region_name=_REGION).create_bucket(
            Bucket=_BUCKET, CreateBucketConfiguration={"LocationConstraint": _REGION}
        )
        yield _BUCKET


def _use_s3(settings, bucket_name):
    settings.MEDIA_STORAGE = "S3"
    settings.STORAGES = {
        **settings.STORAGES,
        "default": {
            "BACKEND": _S3_BACKEND,
            "OPTIONS": {"bucket_name": bucket_name, "region_name": _REGION},
        },
    }


# Settings -------------------------------------------------------------------


def test_filesystem_storage_is_the_default(monkeypatch):
    s = _load_settings(monkeypatch)

    assert s.MEDIA_STORAGE == "FILESYSTEM"
    assert s.STORAGES["default"]["BACKEND"] == _FILESYSTEM_BACKEND
    assert (
        s.STORAGES["staticfiles"]["BACKEND"]
        == "django.contrib.staticfiles.storage.StaticFilesStorage"
    )


def test_s3_storage_takes_its_options_from_the_environment(monkeypatch):
    s = _load_settings(
        monkeypatch,
        MEDIA_STORAGE="s3",
        MEDIA_S3_BUCKET_NAME=_BUCKET,
        MEDIA_S3_LOCATION="media",
        MEDIA_S3_ENDPOINT_URL="https://s3.example.com",
        MEDIA_S3_REGION_NAME=_REGION,
        MEDIA_S3_ACCESS_KEY_ID="key-id",
        MEDIA_S3_SECRET_ACCESS_KEY="secret",
        MEDIA_S3_URL_EXPIRY_S="600",
    )

    assert s.MEDIA_STORAGE == "S3"
    assert s.STORAGES["default"] == {
        "BACKEND": _S3_BACKEND,
        "OPTIONS": {
            "bucket_name": _BUCKET,
            "location": "media",
            "endpoint_url": "https://s3.example.com",
            "region_name": _REGION,
            "access_key_id": "key-id",
            "secret_access_key": "secret",
            "url_expiry": 600,
        },
    }


def test_empty_s3_options_are_left_to_boto3(monkeypatch):
    s = _load_settings(monkeypatch, MEDIA_STORAGE="S3", MEDIA_S3_BUCKET_NAME=_BUCKET)

    options = s.STORAGES["default"]["OPTIONS"]
    assert options["location"] == ""
    assert options["endpoint_url"] is None
    assert options["region_name"] is None
    assert options["access_key_id"] is None
    assert options["secret_access_key"] is None
    assert options["url_expiry"] == 3600


# System check ---------------------------------------------------------------


def test_the_check_is_registered():
    assert check_media_storage in registry.registry.get_checks()


def test_filesystem_storage_passes_the_check(settings):
    settings.MEDIA_STORAGE = "FILESYSTEM"

    assert check_media_storage(None) == []


def test_an_unknown_media_storage_fails_the_check(settings):
    settings.MEDIA_STORAGE = "NFS"

    errors = check_media_storage(None)

    assert [e.id for e in errors] == ["viewer.E001"]
    assert all(isinstance(e, Error) for e in errors)


def test_s3_without_a_bucket_name_fails_the_check(settings):
    _use_s3(settings, "")

    assert [e.id for e in check_media_storage(None)] == ["viewer.E002"]


def test_s3_with_an_existing_bucket_passes_the_check(settings, s3_bucket):
    _use_s3(settings, s3_bucket)

    assert check_media_storage(None) == []


def test_s3_with_a_missing_bucket_fails_the_check(
    settings, s3_bucket
):  # pylint: disable=unused-argument
    _use_s3(settings, "no-such-bucket")

    errors = check_media_storage(None)

    assert [e.id for e in errors] == ["viewer.E003"]
    assert "no-such-bucket" in errors[0].msg
