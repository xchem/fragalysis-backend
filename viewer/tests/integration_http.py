"""Shared HTTP plumbing for the API integration tests.

The integration tests drive a **running** Fragalysis stack purely over HTTP -
they own no Django test DB; the stack owns the one real database. Everything
here is the machinery that is identical whichever endpoint is under test: the
``INTEGRATION_BASE_URL`` gate, a pool manager with sane unattended-CI timeouts,
JSON GETs, paginated result collection and the task-status poll.

Endpoint-specific behaviour (what to POST, what the response looks like, what
to assert) stays in the individual test modules - see
``test_api_upload_target_experiments_async.py`` (LHS) and
``test_viewer_upload_cset_async.py`` (RHS).
"""

import json
import os
import time
from typing import List, Optional

import pytest
import urllib3

#: Base URL of the running stack. Absent on the host/CI run, so the tests skip.
BASE_URL_ENV = "INTEGRATION_BASE_URL"

#: How long to wait for an async task to finish, and how often to poll.
#: Generous because the work runs on a single Celery worker and a slow/contended
#: CI runner has been seen to take ~30 min for the target load alone - 45 min
#: leaves headroom so a sluggish runner does not flake the test on a timeout.
POLL_TIMEOUT_SECONDS = 45 * 60
POLL_INTERVAL_SECONDS = 5

#: Per-request read timeout (seconds). The whole load completes in a few minutes
#: and every GET is sub-second, so a request that goes this long with no data is
#: wedged - fail loudly rather than hang the unattended CI job indefinitely.
REQUEST_READ_TIMEOUT_SECONDS = 120.0

requires_base_url = pytest.mark.skipif(
    not os.environ.get(BASE_URL_ENV),
    reason=f"{BASE_URL_ENV} unset; only the integration stack runs this test.",
)


def base_url() -> str:
    """The running stack's base URL, no trailing slash."""
    return os.environ[BASE_URL_ENV].rstrip("/")


def pool_manager() -> urllib3.PoolManager:
    """A pool manager configured for unattended runs.

    A read timeout is essential: every request is unattended in CI, so a slow
    or wedged endpoint must fail the test loudly rather than hang the job for
    hours. The read timeout is the gap allowed *between* bytes, so it does not
    penalise a large-but-steady upload; retries are off so a stuck request
    surfaces immediately instead of being silently retried.
    """
    return urllib3.PoolManager(
        timeout=urllib3.Timeout(connect=15.0, read=REQUEST_READ_TIMEOUT_SECONDS),
        retries=False,
    )


def get_json(http: urllib3.PoolManager, url: str, **kwargs) -> dict:
    """GET ``url`` and return the decoded JSON body, raising on non-200."""
    response = http.request("GET", url, **kwargs)
    if response.status != 200:
        raise RuntimeError(f"GET {url} returned HTTP {response.status}")
    return json.loads(response.data.decode("utf-8"))


def get_all_results(http: urllib3.PoolManager, url: str, **kwargs) -> list:
    """Return every ``results`` row from a paginated DRF list endpoint.

    Content assertions must search the whole result set, not just the first
    page, so follow the absolute ``next`` link the pagination emits until it is
    null. ``next`` URLs are absolute, so they are requested verbatim.
    """
    results: List[dict] = []
    next_url: Optional[str] = url
    while next_url:
        body = get_json(http, next_url, **kwargs)
        results.extend(body["results"])
        next_url = body.get("next")
    return results


def poll_until_finished(http: urllib3.PoolManager, status_url: str) -> dict:
    """Poll ``status_url`` until the task reports finished, or time out.

    This is the ``api.tasks`` status shape - a ``finished`` flag alongside a
    ``status``. Transient non-200s are expected *while the load runs* and must
    not fail the poll: task_status briefly returns 404 ("Proposal not found")
    because the proposal's Project is not committed until partway through the
    load, and it reports an unfinished body before that. Either just means "not
    ready yet", so keep polling until the deadline; only a finished body ends
    the wait.
    """
    return _poll(
        http,
        status_url,
        is_finished=lambda body: bool(body.get("finished")),
    )


def poll_until_task_state(
    http: urllib3.PoolManager, status_url: str, status_key: str
) -> dict:
    """Poll ``status_url`` until ``status_key`` reaches a terminal Celery state.

    The ``viewer`` task views (``UploadTaskView``, ``ValidateTaskView``) report
    the raw Celery state under their own key (``upload_task_status``,
    ``validate_task_status``) and carry no ``finished`` flag, so "done" means
    the state is terminal - ``SUCCESS`` or ``FAILURE``. The caller decides
    whether the terminal state it got is the one it wanted.
    """
    return _poll(
        http,
        status_url,
        is_finished=lambda body: body.get(status_key) in ("SUCCESS", "FAILURE"),
    )


def _poll(http: urllib3.PoolManager, status_url: str, is_finished) -> dict:
    """Poll ``status_url`` until ``is_finished`` accepts the body, or time out."""
    deadline = time.monotonic() + POLL_TIMEOUT_SECONDS
    last_seen: object = None
    while time.monotonic() < deadline:
        response = http.request("GET", status_url)
        if response.status == 200:
            body = json.loads(response.data.decode("utf-8"))
            if is_finished(body):
                return body
            last_seen = body
        else:
            last_seen = f"HTTP {response.status}"
        time.sleep(POLL_INTERVAL_SECONDS)
    raise AssertionError(
        f"Task at {status_url} did not finish within {POLL_TIMEOUT_SECONDS}s; "
        f"last seen: {last_seen}"
    )
