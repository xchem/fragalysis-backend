"""Fixtures for the frag_s3_storage tests.

S3 is simulated with moto, so these tests need no network and no database.
They deliberately use no fixtures from the backend's apps, so that the module
can later move to its own repository unchanged.
"""
import boto3
import pytest
from moto import mock_aws

from frag_s3_storage import S3Storage

BUCKET = "frag-test-bucket"
REGION = "eu-west-2"


@pytest.fixture
def s3_client(monkeypatch):
    """A moto-backed S3 client with an empty test bucket."""
    # Fake credentials, so that nothing can reach a real AWS account.
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    with mock_aws():
        client = boto3.client("s3", region_name=REGION)
        client.create_bucket(
            Bucket=BUCKET, CreateBucketConfiguration={"LocationConstraint": REGION}
        )
        yield client


@pytest.fixture
def storage(s3_client):  # pylint: disable=redefined-outer-name,unused-argument
    """An S3Storage that uses the moto test bucket."""
    return S3Storage(bucket_name=BUCKET, region_name=REGION)
