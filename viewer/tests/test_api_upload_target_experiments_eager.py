"""Eager (in-process) smoke test for the target-experiment upload API.

Drives a real target archive through the **real** DRF API in the existing
pytest job: download the bundle from the public-read S3 bucket, POST it to
``/api/upload_target_experiments/``, then assert the GET endpoints return the
expected data.

Under pytest the broker runs eagerly (``CELERY_TASK_ALWAYS_EAGER``), so the
view's ``task_load_target.delay(...)`` executes **inline** during the POST and
the database is fully populated by the time the request returns. The eager
backend does not populate the Celery result backend, so this test does *not*
assert ``task_status == SUCCESS`` - that path belongs to the async integration
test. The host needs ``tar`` and ``pigz`` (the loader shells out to them).

Everything here is gated on ``@requires_external_data``; with the env vars
unset the whole module is skipped and the default ``poetry run pytest`` run is
unaffected.
"""
# A test reuses its fixture's name as an argument - a standard pytest pattern
# that pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument
from rest_framework import status

from viewer.tests.external_data import (
    data_identifier,
    download,
    load_manifest,
    requires_external_data,
)

ENDPOINT = "upload_target_experiments"
UPLOAD_URL = "/api/upload_target_experiments/"

#: Maps the manifest `expect` count keys to their GET list endpoints. Counts
#: are only asserted when the manifest value is non-null.
COUNT_ENDPOINTS = {
    "experiments": "/api/experiments/",
    "site_observations": "/api/site_observations/",
    "poses": "/api/poses/",
    "target_experiment_uploads": "/api/target_experiment_uploads/",
}


@requires_external_data
def test_upload_then_get(authenticated_client, settings, tmp_path):
    """Upload the manifest's bundle(s), then assert the GET endpoints match."""
    # No proposal membership exists for the freshly created test user, so turn
    # off upload authorisation (the view returns early when this is False).
    settings.AUTHENTICATE_UPLOAD = False

    entry = load_manifest(ENDPOINT)
    identifier = data_identifier()
    expect = entry["expect"]

    tas_list = []
    for upload in entry["uploads"]:
        tas = upload["tas"]
        filename = upload["file"]
        tas_list.append(tas)

        # Flat S3 layout: <root>/api/<endpoint>/<IDENTIFIER>/<file> - the file
        # sits directly under the identifier, there is no per-TAS directory.
        relative_key = f"api/{ENDPOINT}/{identifier}/{filename}"
        local_file = download(relative_key, tmp_path / filename)

        with open(local_file, "rb") as bundle:
            response = authenticated_client.post(
                UPLOAD_URL,
                data={"target_access_string": tas, "file": bundle},
                format="multipart",
            )

        assert response.status_code == status.HTTP_202_ACCEPTED, response.data
        assert "task_status_url" in response.data

    # The loaded Project.title equals the TAS; publishing the TAS makes the data
    # visible to GET without proposal membership.
    settings.PUBLIC_TAS_LIST = tas_list

    targets = authenticated_client.get("/api/targets/")
    assert targets.status_code == status.HTTP_200_OK
    titles = {row["title"] for row in targets.data["results"]}
    assert expect["target_title"] in titles
    assert targets.data["count"] == expect["targets"]

    for key, url in COUNT_ENDPOINTS.items():
        expected_count = expect.get(key)
        if expected_count is None:
            continue
        response = authenticated_client.get(url)
        assert response.status_code == status.HTTP_200_OK, url
        assert response.data["count"] == expected_count, url
