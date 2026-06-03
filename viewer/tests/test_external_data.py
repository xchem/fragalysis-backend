"""Unit tests for the external-data helper's pure logic.

These run unconditionally (no S3 access, no env vars needed); the network-bound
helpers are exercised by the eager/async tests that the env vars gate on.
"""

import pytest

from viewer.tests import external_data


@pytest.mark.parametrize(
    "s3_url, expected",
    [
        (
            "s3://im-fragalysis-backend/unit-test",
            "https://im-fragalysis-backend.s3.amazonaws.com/unit-test",
        ),
        (
            "s3://im-fragalysis-backend/unit-test/api/x.tgz",
            "https://im-fragalysis-backend.s3.amazonaws.com/unit-test/api/x.tgz",
        ),
        # A bare bucket (no key) yields the bucket root.
        ("s3://bucket-only", "https://bucket-only.s3.amazonaws.com"),
    ],
)
def test_s3_url_to_https(s3_url, expected):
    assert external_data.s3_url_to_https(s3_url) == expected


@pytest.mark.parametrize("bad_url", ["https://example.com/x", "s3://", "not-a-url"])
def test_s3_url_to_https_rejects_non_s3(bad_url):
    with pytest.raises(ValueError):
        external_data.s3_url_to_https(bad_url)


def test_external_data_disabled_without_env(monkeypatch):
    monkeypatch.delenv(external_data.BUCKET_AND_PATH_ENV, raising=False)
    monkeypatch.delenv(external_data.DATA_IDENTIFIER_ENV, raising=False)
    assert external_data.external_data_enabled() is False


def test_external_data_enabled_with_both_env(monkeypatch):
    monkeypatch.setenv(external_data.BUCKET_AND_PATH_ENV, "s3://bucket/prefix")
    monkeypatch.setenv(external_data.DATA_IDENTIFIER_ENV, "ALPHA")
    assert external_data.external_data_enabled() is True


def test_external_data_disabled_with_only_one_env(monkeypatch):
    monkeypatch.setenv(external_data.BUCKET_AND_PATH_ENV, "s3://bucket/prefix")
    monkeypatch.delenv(external_data.DATA_IDENTIFIER_ENV, raising=False)
    assert external_data.external_data_enabled() is False
