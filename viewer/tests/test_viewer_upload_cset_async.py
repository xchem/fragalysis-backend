"""Async integration test for the computed-set (RHS) upload endpoint.

This is the RHS half of the integration story (issue #1014): once the LHS target
archive has loaded (``test_api_upload_target_experiments_async.py``), upload a
computed set against it and prove it processes and becomes publicly visible.

The two tests are **ordered** - this one asserts its precondition (the
manifest's target is present in ``/api/targets/``) up front and fails with a
clear message if the LHS load has not run, because RHS data is meaningless
without it: ``viewer/cset_upload.py`` resolves each molecule's ``ref_mols`` to a
``SiteObservation.code`` produced by that specific load.

``/viewer/upload_cset/`` is a browser-facing view, so this test works rather
differently from the LHS one:

* it returns **HTML**, not JSON. On success the page emits
  ``var taskUrl = "/viewer/upload_task/<uuid>/";`` which has to be scraped out.
* errors are **302 redirects** carrying the message in the session, not error
  statuses - so "no task URL in the page" is the failure signal, and the test
  re-fetches the page to surface the error text.
* ``UploadTaskView`` reports ``upload_task_status``/``upload_task_id``, *not*
  the ``finished``/``status`` pair the LHS poller expects, so it polls to a
  terminal Celery state and then checks ``validated == 'Validated'``.

The upload needs a real logged-in user (the id becomes ``ComputedSet.owner_user``),
so it authenticates as the stack's superuser over HTTP Basic - already in DRF's
``DEFAULT_AUTHENTICATION_CLASSES``, and free of the CSRF dance an admin session
login would need. The GET assertions afterwards are deliberately **anonymous**,
as the LHS test's are: ``PUBLIC_TAS`` publishes the proposal, so this also
proves the computed set is publicly visible.
"""

import os
import re
import warnings

import pytest
import urllib3

from viewer.tests.external_data import (
    download,
    load_manifest,
    missing_objects,
    relative_key,
    requires_external_data,
)
from viewer.tests.integration_http import base_url as _base_url
from viewer.tests.integration_http import get_all_results as _get_all_results
from viewer.tests.integration_http import get_json as _get_json
from viewer.tests.integration_http import (
    poll_until_task_state as _poll_until_task_state,
)
from viewer.tests.integration_http import pool_manager, requires_base_url

#: This test's endpoint, named by its path relative to the test-data root. The
#: same string is the manifest directory and the S3 key prefix.
ENDPOINT = "viewer/upload_cset"

#: Credentials for the upload. The integration stack pins these; the values here
#: are launch-stack.sh's fallback, so a stack that leaves them unset still works.
SUPERUSER_NAME_ENV = "WEB_DJANGO_SUPERUSER_NAME"
SUPERUSER_PASSWORD_ENV = "WEB_DJANGO_SUPERUSER_PASSWORD"

#: The upload page emits the task URL into a script block; this pulls it out.
_RE_TASK_URL = re.compile(r'var taskUrl = "([^"]+)"')

pytestmark = pytest.mark.integration


def _credentials() -> tuple:
    """The superuser username/password the upload authenticates as."""
    return (
        os.environ.get(SUPERUSER_NAME_ENV) or "admin",
        os.environ.get(SUPERUSER_PASSWORD_ENV) or "UNSECURED",
    )


def _assert_lhs_loaded(http: urllib3.PoolManager, base_url: str, target: str) -> None:
    """Fail early and clearly if the LHS load this RHS data depends on is absent."""
    titles = {
        row["title"] for row in _get_all_results(http, f"{base_url}/api/targets/")
    }
    assert target in titles, (
        f"Target {target!r} is not loaded, so the computed set cannot resolve its "
        f"ref_mols. Run the LHS test (test_api_upload_target_experiments_async.py) "
        f"first; /api/targets/ currently holds {sorted(titles)}"
    )


def _post_computed_set(
    http: urllib3.PoolManager,
    base_url: str,
    upload: dict,
    sdf_bytes: bytes,
) -> str:
    """POST the SDF and return the upload task's status URL.

    ``submit_choice=U`` chains validate + process, which is the whole upload in
    one request. ``pdb_zip`` is optional and deliberately not sent - the SDF's
    ``ref_pdb`` values are site-observation codes resolved against the LHS load.
    """
    filename = upload["file"]
    response = http.request(
        "POST",
        f"{base_url}/viewer/upload_cset/",
        headers=urllib3.make_headers(basic_auth=":".join(_credentials())),
        fields={
            "target_name": upload["target"],
            "proposal_ref": upload["tas"],
            "submit_choice": "U",
            "update_set": "None",
            "sdf_file": (filename, sdf_bytes, "chemical/x-mdl-sdfile"),
        },
        # Errors are 302s back to the form; follow them by hand so the redirect
        # itself is visible as the failure signal rather than being swallowed.
        redirect=False,
    )
    assert response.status in (
        200,
        302,
    ), f"upload returned HTTP {response.status}: {response.data[:500]!r}"

    page = response.data.decode("utf-8", "replace")
    match = _RE_TASK_URL.search(page)
    if match:
        return match.group(1)

    # No task URL: the view rejected the request and redirected, stashing the
    # reason in the session. Re-fetch the form to read it back so the assertion
    # names the actual problem rather than "no task URL". urllib3 keeps no
    # cookie jar, so the session cookie has to be carried across by hand -
    # without it the re-fetch gets a *new* session and the message is lost.
    headers = urllib3.make_headers(basic_auth=":".join(_credentials()))
    cookie = response.headers.get("Set-Cookie")
    if cookie:
        headers["Cookie"] = "; ".join(
            part.split(";")[0] for part in cookie.split(", ") if "=" in part
        )
    follow = http.request("GET", f"{base_url}/viewer/upload_cset/", headers=headers)
    reasons = re.findall(
        r"[^<>]*(?:not found|cannot|error)[^<>]*",
        follow.data.decode("utf-8", "replace"),
        re.I,
    )
    raise AssertionError(
        f"upload produced no task URL (HTTP {response.status}, "
        f"Location={response.headers.get('Location')!r}); page said: {reasons[:5]}"
    )


@requires_external_data
@requires_base_url
def test_upload_cset_poll_then_get(tmp_path):
    """Upload a computed set over HTTP, poll to SUCCESS, then assert the GETs."""
    base_url = _base_url()
    http = pool_manager()

    entry = load_manifest(ENDPOINT)
    expect = entry["expect"]

    for upload in entry["uploads"]:
        _assert_lhs_loaded(http, base_url, upload["target"])

        # The LHS site-observation count *before* the upload: each computed
        # molecule adds a virtual SiteObservation, so the growth is a direct
        # check on how many molecules were actually processed.
        before = _get_json(http, f"{base_url}/api/site_observations/")["count"]

        filename = upload["file"]
        local_file = download(relative_key(ENDPOINT, filename), tmp_path / filename)
        with open(local_file, "rb") as sdf:
            status_path = _post_computed_set(http, base_url, upload, sdf.read())

        result = _poll_until_task_state(
            http, f"{base_url}{status_path}", "upload_task_status"
        )
        assert result["upload_task_status"] == "SUCCESS", result.get(
            "upload_traceback", result
        )
        # 'Not validated' carries an HTML table of the validation errors; put it
        # in the message so a failing CI run says *why* the SDF was rejected.
        assert (
            result.get("validated") == "Validated"
        ), f"computed set was not validated: {result.get('html', result)}"

        expected_molecules = expect.get("computed_molecules")
        observed_growth = (
            _get_json(http, f"{base_url}/api/site_observations/")["count"] - before
        )
        if expected_molecules is None:
            warnings.warn(
                f"computed_molecules: manifest count is null; observed "
                f"{observed_growth} new site observations",
                stacklevel=2,
            )
        else:
            assert observed_growth == expected_molecules

    # The GETs are anonymous - PUBLIC_TAS publishes the proposal, so this also
    # proves the computed set is visible without membership.
    compound_sets = _get_json(http, f"{base_url}/api/compound-sets/")
    expected_sets = expect.get("compound_sets")
    if expected_sets is None:
        warnings.warn(
            f"compound_sets: manifest count is null; observed "
            f"{compound_sets['count']} at /api/compound-sets/",
            stacklevel=2,
        )
    else:
        assert compound_sets["count"] == expected_sets

    expected_objects = (expect.get("objects") or {}).get("compound_sets") or []
    if expected_objects:
        results = _get_all_results(http, f"{base_url}/api/compound-sets/")
        missing = missing_objects(results, expected_objects)
        assert (
            not missing
        ), f"expected compound sets not found in /api/compound-sets/: {missing}"
