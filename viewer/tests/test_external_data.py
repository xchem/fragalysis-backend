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


# --- manifest_path / relative_key ------------------------------------------
#
# An endpoint is named by its path RELATIVE to test_data (e.g.
# "api/upload_target_experiments", "viewer/upload_cset"). The same string is
# the manifest directory and the S3 key prefix - that invariant is what lets a
# reader find the bucket objects for a manifest, so it is asserted here.


def test_manifest_path_is_under_test_data():
    path = external_data.manifest_path("api/upload_target_experiments")
    assert path.parent.name == "upload_target_experiments"
    assert path.parent.parent.name == "api"
    assert path.name == "manifest.yaml"


def test_manifest_path_accepts_a_non_api_endpoint():
    path = external_data.manifest_path("viewer/upload_cset")
    assert path.parent.name == "upload_cset"
    assert path.parent.parent.name == "viewer"


def test_relative_key_prefixes_endpoint_and_identifier(monkeypatch):
    monkeypatch.setenv(external_data.DATA_IDENTIFIER_ENV, "ALPHA")
    assert (
        external_data.relative_key("viewer/upload_cset", "compound-set_A71EV2A.sdf")
        == "viewer/upload_cset/ALPHA/compound-set_A71EV2A.sdf"
    )


def test_relative_key_matches_the_manifest_directory(monkeypatch):
    """The S3 key prefix and the manifest directory must be the same string."""
    monkeypatch.setenv(external_data.DATA_IDENTIFIER_ENV, "ALPHA")
    endpoint = "api/upload_target_experiments"
    key = external_data.relative_key(endpoint, "bundle.tgz")
    assert key.startswith(f"{endpoint}/")
    assert external_data.manifest_path(endpoint).parent.match(f"*/{endpoint}")


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
