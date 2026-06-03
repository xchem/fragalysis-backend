"""Async integration test for the target-experiment upload API.

The faithful counterpart to the eager smoke test: instead of running in-process
under an eager broker, this drives a **running** Fragalysis stack (a real DB, a
redis broker and a celery worker) purely over HTTP. It uploads the bundle,
polls the ``task_status`` endpoint until the load finishes, asserts the task
reached ``SUCCESS``, then asserts the GET endpoints return the expected data.

It is marked ``integration`` (deselected by the default ``-m "not integration"``
addopts) **and** gated on ``INTEGRATION_BASE_URL`` (the base URL of the running
stack, e.g. ``http://backend`` - nginx listens on :80 inside the network), so
it is never collected by the host or
CI pytest run - only ``docker-compose.integration.yml`` sets that variable.

Visibility for the anonymous GET/poll comes from the stack's ``PUBLIC_TAS``,
which must include the manifest's TAS; the upload itself runs with the stack's
``AUTHENTICATE_UPLOAD`` off.
"""
import json
import os
import time

import pytest
import urllib3

from viewer.tests.external_data import (
    data_identifier,
    download,
    load_manifest,
    requires_external_data,
)

ENDPOINT = "upload_target_experiments"

#: Base URL of the running stack. Absent on the host/CI run, so the test skips.
BASE_URL_ENV = "INTEGRATION_BASE_URL"

#: How long to wait for the async load to finish, and how often to poll.
POLL_TIMEOUT_SECONDS = 30 * 60
POLL_INTERVAL_SECONDS = 5

COUNT_ENDPOINTS = {
    "experiments": "/api/experiments/",
    "site_observations": "/api/site_observations/",
    "poses": "/api/poses/",
    "target_experiment_uploads": "/api/target_experiment_uploads/",
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


def _poll_until_finished(http: urllib3.PoolManager, status_url: str) -> dict:
    """Poll ``status_url`` until the task reports finished, or time out."""
    deadline = time.monotonic() + POLL_TIMEOUT_SECONDS
    while True:
        body = _get_json(http, status_url)
        if body.get("finished"):
            return body
        if time.monotonic() >= deadline:
            raise AssertionError(
                f"Task at {status_url} did not finish within "
                f"{POLL_TIMEOUT_SECONDS}s; last status: {body}"
            )
        time.sleep(POLL_INTERVAL_SECONDS)


@requires_external_data
@requires_base_url
def test_upload_poll_then_get(tmp_path):
    """Upload over HTTP, poll to SUCCESS, then assert the GET endpoints match."""
    base_url = _base_url()
    http = urllib3.PoolManager()

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
        expected_count = expect.get(key)
        if expected_count is None:
            continue
        body = _get_json(http, f"{base_url}{url}")
        assert body["count"] == expected_count, url
