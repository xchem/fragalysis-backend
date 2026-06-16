"""Async integration test for the target-experiment upload API.

The target-upload flow is genuinely asynchronous (the API hands the load to a
Celery worker), so it is only exercised faithfully against a **running**
Fragalysis stack - a real DB, a redis broker and a celery worker - driven purely
over HTTP. Running it under an eager, in-process broker is unwise: eager mode
does not model the broker/worker handshake the upload view performs, so this
test deliberately does *not* have an in-process counterpart. It uploads the
bundle, polls the ``task_status`` endpoint until the load finishes, asserts the
task reached ``SUCCESS``, then asserts the GET endpoints return the expected
data. Finally - per issue #967 - it builds and fetches a download archive for
the loaded target, proving the download endpoint produces a real, non-empty zip
(the archive itself is not kept).

It is marked ``integration`` (deselected by the default ``-m "not integration"``
addopts) **and** gated on ``INTEGRATION_BASE_URL`` (the base URL of the running
stack, e.g. ``http://backend`` - nginx listens on :80 inside the network), so
it is never collected by the host or CI pytest run - only
``docker-compose.integration.yml`` sets that variable.

Visibility for the anonymous GET/poll comes from the stack's ``PUBLIC_TAS``,
which must include the manifest's TAS; the upload itself runs with the stack's
``AUTHENTICATE_UPLOAD`` off.
"""
import json
import os
import time
import warnings
from typing import List, Optional

import pytest
import urllib3

from viewer.tests.external_data import (
    data_identifier,
    download,
    load_manifest,
    missing_objects,
    requires_external_data,
)

ENDPOINT = "upload_target_experiments"

#: Base URL of the running stack. Absent on the host/CI run, so the test skips.
BASE_URL_ENV = "INTEGRATION_BASE_URL"

#: How long to wait for the async load to finish, and how often to poll.
POLL_TIMEOUT_SECONDS = 30 * 60
POLL_INTERVAL_SECONDS = 5

#: Per-request read timeout (seconds). The whole load completes in a few minutes
#: and every GET is sub-second, so a request that goes this long with no data is
#: wedged - fail loudly rather than hang the unattended CI job indefinitely.
REQUEST_READ_TIMEOUT_SECONDS = 120.0

#: Endpoints whose ``count`` is asserted against the manifest. These are the
#: model objects a target load produces, beyond the target itself. Each key
#: must match a manifest ``expect`` key; a null/absent value turns the
#: assertion off (see the manifest header), so new endpoints can land here and
#: gain assertions as the real counts are filled in.
COUNT_ENDPOINTS = {
    "experiments": "/api/experiments/",
    "site_observations": "/api/site_observations/",
    "canon_sites": "/api/canon_sites/",
    "canon_site_confs": "/api/canon_site_confs/",
    "xtalform_sites": "/api/xtalform_sites/",
    "poses": "/api/poses/",
    "target_experiment_uploads": "/api/target_experiment_uploads/",
}

#: Endpoints whose returned *objects* can be content-checked against the
#: manifest ``expect.objects`` section (each entry is a field-subset that must
#: match some returned object). Keyed by the manifest endpoint name.
OBJECT_ENDPOINTS = {
    "targets": "/api/targets/",
    "experiments": "/api/experiments/",
    "site_observations": "/api/site_observations/",
    "canon_sites": "/api/canon_sites/",
    "canon_site_confs": "/api/canon_site_confs/",
    "xtalform_sites": "/api/xtalform_sites/",
    "poses": "/api/poses/",
}

pytestmark = pytest.mark.integration

requires_base_url = pytest.mark.skipif(
    not os.environ.get(BASE_URL_ENV),
    reason=f"{BASE_URL_ENV} unset; only the integration stack runs this test.",
)


def _base_url() -> str:
    return os.environ[BASE_URL_ENV].rstrip("/")


def _get_json(http: urllib3.PoolManager, url: str) -> dict:
    """GET ``url`` and return the decoded JSON body, raising on non-200."""
    response = http.request("GET", url)
    if response.status != 200:
        raise RuntimeError(f"GET {url} returned HTTP {response.status}")
    return json.loads(response.data.decode("utf-8"))


def _get_all_results(http: urllib3.PoolManager, url: str) -> list:
    """Return every ``results`` row from a paginated DRF list endpoint.

    Content assertions must search the whole result set, not just the first
    page, so follow the absolute ``next`` link the pagination emits until it is
    null. ``next`` URLs are absolute, so they are requested verbatim.
    """
    results: List[dict] = []
    next_url: Optional[str] = url
    while next_url:
        body = _get_json(http, next_url)
        results.extend(body["results"])
        next_url = body.get("next")
    return results


def _poll_until_finished(http: urllib3.PoolManager, status_url: str) -> dict:
    """Poll ``status_url`` until the task reports finished, or time out.

    Transient non-200s are expected *while the load runs* and must not fail the
    poll: task_status briefly returns 404 ("Proposal not found") because the
    proposal's Project is not committed until partway through the load, and it
    reports an unfinished body before that. Either just means "not ready yet",
    so keep polling until the deadline; only a finished body ends the wait.
    """
    deadline = time.monotonic() + POLL_TIMEOUT_SECONDS
    last_seen: object = None
    while time.monotonic() < deadline:
        response = http.request("GET", status_url)
        if response.status == 200:
            body = json.loads(response.data.decode("utf-8"))
            if body.get("finished"):
                return body
            last_seen = body
        else:
            last_seen = f"HTTP {response.status}"
        time.sleep(POLL_INTERVAL_SECONDS)
    raise AssertionError(
        f"Task at {status_url} did not finish within {POLL_TIMEOUT_SECONDS}s; "
        f"last seen: {last_seen}"
    )


def _download_target(
    http: urllib3.PoolManager, base_url: str, target_name: str, tas: str
) -> None:
    """Build and fetch a download archive for an already-loaded target (#967).

    Mirrors the upload flow: POST to start the build, poll to SUCCESS, then
    re-POST to retrieve the now-ready ``file_url``, then GET the archive itself.
    The file content is not kept - we only assert that a real, non-empty zip
    comes back, which is the whole point of issue #967.
    """
    fields = {"target_name": target_name, "target_access_string": tas}

    response = http.request(
        "POST", f"{base_url}/api/download_structures/", fields=fields
    )
    # 202: a build task was started (or one is already in progress). 200: a
    # matching completed download already exists and the file_url is returned
    # immediately. Both are valid starting points.
    assert response.status in (202, 200), (response.status, response.data)
    body = json.loads(response.data.decode("utf-8"))

    if response.status == 202:
        status_path = body["task_status_url"]
        result = _poll_until_finished(http, f"{base_url}{status_path}")
        assert result["status"] == "SUCCESS", result
        # Re-POST the identical request: the completed DownloadLinks record now
        # exists, so the view returns the ready file_url (200) rather than
        # starting another build. This is more robust than scraping the
        # file_url out of the task status messages.
        response = http.request(
            "POST", f"{base_url}/api/download_structures/", fields=fields
        )
        assert response.status == 200, (response.status, response.data)
        body = json.loads(response.data.decode("utf-8"))

    file_url = body["file_url"]

    archive = http.request(
        "GET", f"{base_url}/api/download_structures/?file_url={file_url}"
    )
    assert archive.status == 200, (archive.status, archive.data[:200])
    # The endpoint labels every download "application/zip", but the actual body
    # depends on the requested format: the default (use_zip off) is a gzip
    # tarball (pigz), and use_zip would yield a real zip. We only need to prove
    # a real, non-empty archive came back, so accept either signature rather
    # than couple this test to the compression choice.
    assert int(archive.headers.get("Content-Length", "0")) > 0
    zip_magic, gzip_magic = b"PK\x03\x04", b"\x1f\x8b"
    assert (
        archive.data[:4] == zip_magic or archive.data[:2] == gzip_magic
    ), f"download body is not a zip/gzip archive: {archive.data[:4]!r}"


@requires_external_data
@requires_base_url
def test_upload_poll_then_get(tmp_path):
    """Upload over HTTP, poll to SUCCESS, then assert the GET endpoints match."""
    base_url = _base_url()
    # A read timeout is essential: every request below is unattended in CI, so a
    # slow or wedged endpoint must fail the test loudly rather than hang the job
    # for hours. The read timeout is the gap allowed *between* bytes, so it does
    # not penalise a large-but-steady upload; retries are off so a stuck request
    # surfaces immediately instead of being silently retried.
    http = urllib3.PoolManager(
        timeout=urllib3.Timeout(connect=15.0, read=REQUEST_READ_TIMEOUT_SECONDS),
        retries=False,
    )

    entry = load_manifest(ENDPOINT)
    identifier = data_identifier()
    expect = entry["expect"]

    for upload in entry["uploads"]:
        tas = upload["tas"]
        filename = upload["file"]

        relative_key = f"api/{ENDPOINT}/{identifier}/{filename}"
        local_file = download(relative_key, tmp_path / filename)

        with open(local_file, "rb") as bundle:
            response = http.request(
                "POST",
                f"{base_url}/api/upload_target_experiments/",
                fields={
                    "target_access_string": tas,
                    "file": (filename, bundle.read(), "application/gzip"),
                },
            )
        assert (
            response.status == 202
        ), f"upload returned HTTP {response.status}: {response.data!r}"

        status_path = json.loads(response.data.decode("utf-8"))["task_status_url"]
        result = _poll_until_finished(http, f"{base_url}{status_path}")
        assert result["status"] == "SUCCESS", result

    targets = _get_json(http, f"{base_url}/api/targets/")
    titles = {row["title"] for row in targets["results"]}
    assert expect["target_title"] in titles
    assert targets["count"] == expect["targets"]

    for key, url in COUNT_ENDPOINTS.items():
        body = _get_json(http, f"{base_url}{url}")
        observed = body["count"]
        expected_count = expect.get(key)
        if expected_count is None:
            # The manifest has no count for this endpoint yet. Don't assert -
            # surface the observed count so it can be filled into the manifest
            # (see "Adding a new dataset" in INTEGRATION-TESTS.md). The warning
            # shows in pytest's summary even on an otherwise-passing run, so a
            # single CI run reveals every count waiting to be pinned down.
            warnings.warn(
                f"{key}: manifest count is null; observed {observed} at {url}",
                stacklevel=2,
            )
            continue
        assert observed == expected_count, url

    # Content assertions: beyond the bare counts, the manifest can name specific
    # objects (as field-subsets) it expects an endpoint to return for this
    # identifier. Scan every page and fail loudly with the unmatched
    # expectations if any are missing.
    expected_objects = expect.get("objects") or {}
    for key, expected_list in expected_objects.items():
        if not expected_list:
            continue
        object_url = OBJECT_ENDPOINTS.get(key)
        assert (
            object_url is not None
        ), f"manifest 'objects' names unknown endpoint {key!r}"
        results = _get_all_results(http, f"{base_url}{object_url}")
        missing = missing_objects(results, expected_list)
        assert (
            not missing
        ), f"{key}: expected objects not found in {object_url}: {missing}"

    # Finally (#967): prove the loaded target can be downloaded. Reuse the
    # already-loaded state rather than re-uploading - the upload above is the
    # expensive part. The manifest already carries the title and TAS we need.
    _download_target(http, base_url, expect["target_title"], entry["uploads"][0]["tas"])
