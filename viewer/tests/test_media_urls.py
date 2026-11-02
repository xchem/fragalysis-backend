# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the file URLs that the API returns (issue 1061).

The browser must get ``/media/...`` URLs in both media storages. Django checks
access for these URLs, then nginx serves the file (see ``viewer.media_storage``).
In S3 mode, a presigned S3 URL is for nginx only. The API must never return one,
because it is cross-origin, it expires, and anyone who has it can use it.
"""
import pytest
from django.test import RequestFactory
from rest_framework import serializers

from viewer.models import Target

_TARGET_DIR = "Mpro_lb1234"
_METADATA = "metadata/Mpro meta#1.csv"
_METADATA_URL_PATH = "/media/metadata/Mpro%20meta%231.csv"


@pytest.fixture
def member_target(media_backend, user, make_project, make_target, settings):
    settings.TA_AUTH_SERVICE = ""
    target = make_target(make_project("proposal-urls", members=[user]))
    target.zip_archive = _TARGET_DIR
    target.metadata = _METADATA
    target.save()
    return target


def _only_target(response):
    assert response.status_code == 200
    (target,) = response.json()["results"]
    return target


def _assert_not_s3(value):
    assert "X-Amz" not in value
    assert "Signature" not in value
    assert "frag-media-test" not in value


def test_targets_api_returns_the_target_directory_name(
    authenticated_client, member_target
):
    response = authenticated_client.get("/api/targets/")

    assert _only_target(response)["zip_archive"] == _TARGET_DIR


def test_targets_api_returns_null_for_no_target_directory(
    authenticated_client, member_target
):
    member_target.zip_archive = None
    member_target.save()

    response = authenticated_client.get("/api/targets/")

    assert _only_target(response)["zip_archive"] is None


def test_targets_api_returns_a_media_url_for_the_metadata(
    authenticated_client, member_target
):
    response = authenticated_client.get("/api/targets/")

    metadata = _only_target(response)["metadata"]
    assert metadata == "https://testserver" + _METADATA_URL_PATH
    _assert_not_s3(metadata)


class _TargetFilesSerializer(serializers.ModelSerializer):
    """A plain DRF serializer, like the ``fields = "__all__"`` ones in the API."""

    class Meta:
        model = Target
        fields = ["metadata"]


def test_drf_file_fields_are_media_urls(member_target):
    request = RequestFactory().get("/api/")

    data = _TargetFilesSerializer(member_target, context={"request": request}).data

    assert data["metadata"] == "http://testserver" + _METADATA_URL_PATH
    _assert_not_s3(data["metadata"])
