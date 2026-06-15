"""Broad smoke tests of the REST API against an *empty* installation.

The database here has no ``Target`` (and therefore none of the experiment /
site-observation / pose data that hangs off one). The point is breadth, not
depth: confirm that the API as a whole behaves sanely with nothing loaded -
list endpoints answer 200 with an empty, paginated result set; access control
still applies; read-only viewsets stay read-only; and the supporting
user/role endpoints respond as documented. Deep, data-driven assertions belong
with the fixtures that load a real target archive, not here.
"""

import pytest

# A deliberately broad slice of the access-controlled (ISPyBSafeQuerySet) list
# endpoints. With an empty database each must answer 200 with an empty page,
# whether or not the caller is authenticated - there is simply nothing to show.
LIST_ENDPOINTS = [
    "/api/targets/",
    "/api/projects/",
    "/api/compounds/",
    "/api/experiments/",
    "/api/target_experiment_uploads/",
    "/api/site_observations/",
    "/api/canon_sites/",
    "/api/canon_site_confs/",
    "/api/xtalform_sites/",
    "/api/poses/",
    "/api/session-projects/",
    "/api/snapshots/",
    "/api/compound-sets/",
    "/api/tag_category/",
    "/api/siteobservation_tag/",
    "/api/session_project_tag/",
]


def _results(response):
    """The result list from a LimitOffset-paginated response."""
    assert "results" in response.data, response.data
    return response.data["results"]


@pytest.mark.parametrize("endpoint", LIST_ENDPOINTS)
def test_list_endpoint_empty_for_authenticated_user(authenticated_client, endpoint):
    """Every list endpoint is reachable and empty for an authenticated user."""
    response = authenticated_client.get(endpoint)

    assert response.status_code == 200
    assert _results(response) == []
    assert response.data["count"] == 0


@pytest.mark.django_db
@pytest.mark.parametrize("endpoint", LIST_ENDPOINTS)
def test_list_endpoint_empty_for_anonymous_user(api_client, endpoint):
    """The same endpoints are reachable and empty for an anonymous caller.

    Access control yields only public data to an anonymous user; with nothing
    loaded that is, correctly, an empty page rather than an error.
    """
    response = api_client.get(endpoint)

    assert response.status_code == 200
    assert _results(response) == []


def test_retrieve_missing_target_is_404(authenticated_client):
    """Retrieving a non-existent record returns 404, not a server error."""
    response = authenticated_client.get("/api/targets/999999/")

    assert response.status_code == 404


def test_readonly_viewset_rejects_post(authenticated_client):
    """Projects are read-only over the API: POST is not allowed (405)."""
    response = authenticated_client.post("/api/projects/", data={"title": "nope"})

    assert response.status_code == 405


def test_user_endpoint_reports_no_access_for_authenticated_user(authenticated_client):
    """``/api/user/`` describes the caller and, with no service, no access.

    With ``TA_AUTH_SERVICE`` unset the authenticator reports
    ``SERVICE_NOT_PRESENT`` and the user has no target access - the expected
    shape for an empty, service-less installation.
    """
    response = authenticated_client.get("/api/user/")

    assert response.status_code == 200
    body = response.json()
    assert body["user"] == "tester"
    assert body["target_access"] == []
    assert body["ping"] == "SERVICE_NOT_PRESENT"
    assert body["authenticator"]["kind"] == "SERVICE_NOT_PRESENT"


@pytest.mark.django_db
def test_user_roles_requires_authentication(api_client):
    """``/api/user_roles/`` is authentication-gated."""
    response = api_client.get("/api/user_roles/")

    assert response.status_code in (401, 403)


def test_user_roles_lists_seeded_roles_for_authenticated_user(authenticated_client):
    """The role list is data-independent: the migration-seeded 'Loader' role is
    present even on an empty installation, and holds no users."""
    response = authenticated_client.get("/api/user_roles/")

    assert response.status_code == 200
    roles = {row["name"]: row for row in _results(response)}
    assert "Loader" in roles
    assert roles["Loader"]["user_count"] == 0
