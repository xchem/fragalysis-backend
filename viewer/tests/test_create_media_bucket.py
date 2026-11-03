# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the create_media_bucket command (issue 1038).

Compose stacks (development and the integration tests) run an S3-compatible
store (SeaweedFS) that starts empty. A one-shot service runs this command to
create the media bucket before the backend starts, because the start-up check
(viewer.E003) needs the bucket. It works with any S3 endpoint.
"""
from io import StringIO

import boto3
import pytest
from botocore.exceptions import EndpointConnectionError
from django.core.management import CommandError, call_command
from moto import mock_aws

_BUCKET = "frag-media-init"
_REGION = "eu-west-2"


@pytest.fixture
def s3(settings, monkeypatch):
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    with mock_aws():
        settings.MEDIA_S3_STORAGE_OPTIONS = {
            "bucket_name": _BUCKET,
            "region_name": _REGION,
        }
        yield boto3.client("s3", region_name=_REGION)


def _run(*argv) -> str:
    out = StringIO()
    call_command("create_media_bucket", *argv, stdout=out, stderr=out)
    return out.getvalue()


def _buckets(client) -> list[str]:
    return [b["Name"] for b in client.list_buckets()["Buckets"]]


def test_the_bucket_is_created(s3):
    output = _run()

    assert _buckets(s3) == [_BUCKET]
    assert "Created" in output


def test_an_existing_bucket_is_left_alone(s3):
    s3.create_bucket(
        Bucket=_BUCKET, CreateBucketConfiguration={"LocationConstraint": _REGION}
    )

    output = _run()

    assert _buckets(s3) == [_BUCKET]
    assert "already exists" in output


def _mock_client(mocker, head_bucket_side_effect):
    client = mocker.Mock()
    client.head_bucket.side_effect = head_bucket_side_effect
    mocker.patch(
        "frag_s3_storage.storage.S3Storage.client",
        new_callable=mocker.PropertyMock,
        return_value=client,
    )
    return client


def test_the_endpoint_is_waited_for(s3, mocker):
    """The S3 container can still be starting: try again, then succeed."""
    mocker.patch("viewer.management.commands.create_media_bucket.time.sleep")
    client = _mock_client(
        mocker, [EndpointConnectionError(endpoint_url="http://s3:8333"), {}]
    )

    output = _run("--wait-s", "5")

    assert client.head_bucket.call_count == 2
    assert "already exists" in output


def test_an_endpoint_that_never_answers_fails(s3, mocker):
    mocker.patch("viewer.management.commands.create_media_bucket.time.sleep")

    _mock_client(mocker, EndpointConnectionError(endpoint_url="http://s3:8333"))

    with pytest.raises(CommandError, match="could not be reached"):
        _run("--wait-s", "3")


def test_no_bucket_name_fails(settings):
    settings.MEDIA_S3_STORAGE_OPTIONS = {"bucket_name": ""}

    with pytest.raises(CommandError, match="MEDIA_S3_BUCKET_NAME"):
        _run()
