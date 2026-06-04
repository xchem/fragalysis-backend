"""Unit tests for the external-data helper's pure logic.

These run unconditionally (no S3 access, no env vars needed); the network-bound
helpers are exercised by the async integration test that the env vars gate on.
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


# --- object_matches / missing_objects --------------------------------------
#
# These back the integration test's content assertions: a manifest names the
# objects (as field-subsets) it expects each endpoint to return, and these
# helpers decide whether a returned object satisfies an expectation.


def test_object_matches_subset_with_extra_fields():
    # The returned object carries extra fields; only the expected subset matters.
    row = {"code": "x0123a", "compound_code": "Z42", "id": 7}
    assert external_data.object_matches(row, {"code": "x0123a"}) is True
    assert (
        external_data.object_matches(row, {"code": "x0123a", "compound_code": "Z42"})
        is True
    )


def test_object_matches_empty_expectation_matches_anything():
    assert external_data.object_matches({"code": "x0123a"}, {}) is True


def test_object_matches_rejects_missing_key():
    assert (
        external_data.object_matches({"code": "x0123a"}, {"compound_code": "Z42"})
        is False
    )


def test_object_matches_rejects_value_mismatch():
    assert external_data.object_matches({"code": "x0123a"}, {"code": "other"}) is False


def test_missing_objects_returns_empty_when_all_found():
    results = [{"code": "a", "id": 1}, {"code": "b", "id": 2}]
    expected = [{"code": "a"}, {"code": "b"}]
    assert external_data.missing_objects(results, expected) == []


def test_missing_objects_returns_the_unmatched_expectations():
    results = [{"code": "a", "id": 1}]
    expected = [{"code": "a"}, {"code": "b"}, {"code": "c"}]
    assert external_data.missing_objects(results, expected) == [
        {"code": "b"},
        {"code": "c"},
    ]
