"""Shared helpers for tests that pull large data from a public-read S3 bucket.

The API data-graph tests need realistic target archives that are far too large
to commit to the repository. Instead they are stored in a **public-read** S3
bucket (anonymous ``GET`` only; no ``ListBucket``) and downloaded over plain
HTTPS at test time. To keep the default ``poetry run pytest`` run unaffected,
every test that uses this mechanism is gated on two environment variables:

``UNIT_TEST_BUCKET_AND_PATH``
    The ``s3://<bucket>/<prefix>`` root the test data lives under, e.g.
    ``s3://im-fragalysis-backend/unit-test``.

``UNIT_TEST_DATA_IDENTIFIER``
    Selects one entry from a per-endpoint ``manifest.yaml`` (e.g. ``ALPHA``),
    naming the object(s) to upload and the results to expect.

When either variable is absent the ``requires_external_data`` marker skips the
test, mirroring ``requires_archive`` in ``test_target_loader.py``.
"""

import os
from pathlib import Path
from typing import Any, Dict, List
from urllib.parse import urlparse

import pytest
import urllib3
import yaml

#: Environment variable naming the ``s3://bucket/prefix`` test-data root.
BUCKET_AND_PATH_ENV = "UNIT_TEST_BUCKET_AND_PATH"

#: Environment variable selecting the manifest entry (e.g. "ALPHA").
DATA_IDENTIFIER_ENV = "UNIT_TEST_DATA_IDENTIFIER"

#: Where the per-endpoint manifests live, relative to this file.
_TEST_DATA_ROOT = Path(__file__).parent / "test_data" / "api"


def external_data_enabled() -> bool:
    """True when both env vars that enable external-data tests are set."""
    return bool(os.environ.get(BUCKET_AND_PATH_ENV)) and bool(
        os.environ.get(DATA_IDENTIFIER_ENV)
    )


#: Skip marker for tests that need the external S3 test data. Use as a
#: decorator (``@requires_external_data``), exactly like ``requires_archive``.
requires_external_data = pytest.mark.skipif(
    not external_data_enabled(),
    reason=(
        f"External test data disabled; set {BUCKET_AND_PATH_ENV} and "
        f"{DATA_IDENTIFIER_ENV} to enable (see viewer/tests/external_data.py)."
    ),
)


def data_identifier() -> str:
    """The selected manifest entry id (the ``UNIT_TEST_DATA_IDENTIFIER`` value)."""
    return os.environ[DATA_IDENTIFIER_ENV]


def s3_url_to_https(s3_url: str) -> str:
    """Convert an ``s3://bucket/key`` URL to its anonymous HTTPS form.

    e.g. ``s3://im-fragalysis-backend/unit-test`` becomes
    ``https://im-fragalysis-backend.s3.amazonaws.com/unit-test``. The key part
    is optional; a bare ``s3://bucket`` yields the bucket root.
    """
    parsed = urlparse(s3_url)
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError(f"Not an s3:// URL: {s3_url!r}")
    bucket = parsed.netloc
    key = parsed.path.lstrip("/")
    base = f"https://{bucket}.s3.amazonaws.com"
    return f"{base}/{key}" if key else base


def _root_https_url() -> str:
    """The HTTPS form of the ``UNIT_TEST_BUCKET_AND_PATH`` root, no trailing slash."""
    return s3_url_to_https(os.environ[BUCKET_AND_PATH_ENV]).rstrip("/")


def download(relative_key: str, dest: Path) -> Path:
    """Stream an object from ``<root>/<relative_key>`` to ``dest`` (anonymous GET).

    Raises ``RuntimeError`` on any non-200 response - a missing or unreadable
    object must fail the test loudly rather than silently producing an empty or
    partial file.
    """
    url = f"{_root_https_url()}/{relative_key.lstrip('/')}"
    http = urllib3.PoolManager()
    with http.request("GET", url, preload_content=False) as response:
        if response.status != 200:
            response.release_conn()
            raise RuntimeError(
                f"Anonymous GET of {url} returned HTTP {response.status} "
                f"(expected 200)"
            )
        with open(dest, "wb") as out_file:
            for chunk in response.stream(1024 * 1024):
                out_file.write(chunk)
        response.release_conn()
    return dest


def load_manifest(endpoint: str) -> Dict[str, Any]:
    """Return the manifest entry for the current ``UNIT_TEST_DATA_IDENTIFIER``.

    ``endpoint`` names the sub-directory under ``test_data/api/`` holding the
    ``manifest.yaml`` (e.g. ``upload_target_experiments``). If the manifest has
    no entry for the selected identifier the test is skipped rather than
    failing, so a bucket carrying only some identifiers stays usable.
    """
    manifest_path = _TEST_DATA_ROOT / endpoint / "manifest.yaml"
    if not manifest_path.is_file():
        pytest.skip(f"No manifest at {manifest_path}")

    with open(manifest_path, "rt", encoding="utf-8") as manifest_file:
        manifest = yaml.safe_load(manifest_file) or {}

    identifier = data_identifier()
    if identifier not in manifest:
        pytest.skip(
            f"Manifest {manifest_path} has no entry for "
            f"{DATA_IDENTIFIER_ENV}={identifier!r}"
        )
    return manifest[identifier]


def object_matches(row: Dict[str, Any], expected: Dict[str, Any]) -> bool:
    """True when ``row`` satisfies every field in the ``expected`` subset.

    A returned API object carries far more fields than a manifest wants to pin
    down, so an expectation is a *subset*: every key in ``expected`` must be
    present in ``row`` with an equal value, but ``row`` may carry extra fields.
    An empty ``expected`` matches any row.
    """
    return all(key in row and row[key] == value for key, value in expected.items())


def missing_objects(
    results: List[Dict[str, Any]], expected_objects: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Return the ``expected_objects`` that no row in ``results`` satisfies.

    Each expected object is a field-subset (see :func:`object_matches`). An
    empty return value means every expectation was found; a non-empty one names
    exactly which expectations were not, so the caller can fail loudly with the
    unmatched expectations rather than a bare boolean.
    """
    return [
        expected
        for expected in expected_objects
        if not any(object_matches(row, expected) for row in results)
    ]
