"""Tests for the ``/api/tas/`` endpoint (issue #1017).

The endpoint answers "who is a member of this target access string?", by way of
the authenticator's ``/users/{tas}`` endpoint (TA-Auth 1.5.0). It is shaped like
``/api/user/``, which answers the mirror-image question, and carries the same
authenticator/ping block. Two things are worth stating up front, because they
shape every test here:

- **Authentication is the only gate.** Any logged-in user may ask about any TAS,
  whether or not they are a member of it. There is no membership check - the
  tests below assert that a caller with no relationship to the TAS is answered,
  not refused. Only an anonymous caller is turned away.
- **"Nobody" is not "we do not know".** The authenticator answers 503, rather
  than an empty set, when it cannot reach ISPyB, and the client preserves that
  distinction. The endpoint must too - reporting an empty membership for a
  service outage would be a silent lie about who has access.

The authenticator call itself is mocked via ``mock_tas_users``, so no test
contacts a real service.
"""
# The `no_ta_service` fixture is requested for its side effect only - the
# standard pytest pattern that pylint misreads (as it does in conftest.py).
# pylint: disable=unused-argument

import pytest

_TAS = "lb12345-1"


@pytest.fixture(name="no_ta_service")
def fixture_no_ta_service(settings):
    """Keep the TA service out of the picture (it plays no part here)."""
    settings.TA_AUTH_SERVICE = ""


@pytest.mark.django_db
def test_requires_authentication(api_client, mock_tas_users, no_ta_service):
    """An anonymous caller is turned away - being logged in is the one gate."""
    mock_tas_users(users=["abc12345"])

    response = api_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code in (401, 403)


def test_authenticated_user_sees_the_users_of_the_tas(
    authenticated_client, mock_tas_users, no_ta_service
):
    """A logged-in caller gets the membership the authenticator reports."""
    mock_tas_users(users=["xyz98765", "abc12345"])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    body = response.json()
    assert body["tas"] == _TAS
    # Sorted, so the response is stable for the caller (the source is a set).
    assert body["users"] == ["abc12345", "xyz98765"]
    # The whole contract, so the shape cannot drift away from /api/user/'s.
    assert sorted(body.keys()) == ["authenticator", "ping", "tas", "users"]


def test_non_member_may_query_any_tas(
    authenticated_client, make_user, make_project, mock_tas_users, no_ta_service
):
    """Membership is deliberately *not* required to ask.

    The caller here has no relationship with the TAS at all - it belongs to
    somebody else - and is still answered. This is the endpoint's access rule,
    so it is pinned rather than left implicit.
    """
    make_project(_TAS, members=[make_user("someone-else")])
    mock_tas_users(users=["someone-else"])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    assert response.json()["users"] == ["someone-else"]


def test_tas_unknown_to_fragalysis_is_still_queried(
    authenticated_client, mock_tas_users, no_ta_service
):
    """The TAS need not correspond to a Fragalysis Project.

    The question is put to the authenticator (and so to ISPyB), which knows
    about proposals this deployment may never have loaded a target for.
    """
    mock_tas_users(users=["abc12345"])

    response = authenticated_client.get("/api/tas/lb99999-9/")

    assert response.status_code == 200
    body = response.json()
    assert body["tas"] == "lb99999-9"
    assert body["users"] == ["abc12345"]


def test_response_names_the_authenticator_that_answered(
    authenticated_client, mock_tas_users, mock_authenticator, no_ta_service
):
    """The response carries the same authenticator/ping block as /api/user/.

    The two endpoints are a pair - both are answers *from the authenticator* -
    so a caller can tell which service, at which version, produced the list.
    """
    mock_tas_users(users=["abc12345"])
    mock_authenticator(
        version="1.5.0", kind="ISPYB", location="https://ta-auth.example.ac.uk"
    )

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    body = response.json()
    assert body["ping"] == "OK"
    assert body["authenticator"] == {
        "kind": "ISPYB",
        "name": "XChem Python FastAPI TAS Authenticator",
        "version": "1.5.0",
        "location": "https://ta-auth.example.ac.uk",
    }


def test_tas_with_no_members_is_an_empty_list_not_an_error(
    authenticated_client, mock_tas_users, no_ta_service
):
    """An empty membership is a legitimate answer, and answered as 200."""
    mock_tas_users(users=[])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    body = response.json()
    assert body["tas"] == _TAS
    assert body["users"] == []


def test_authenticator_failure_is_reported_not_hidden(
    authenticated_client, mock_tas_users, no_ta_service
):
    """When the authenticator cannot answer, say so - do not report "nobody"."""
    mock_tas_users(error="Status was 503 (not 200)")

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 503
    body = response.json()
    assert body["tas"] == _TAS
    assert "error" in body
    # No 'users' key at all - an empty list here would read as "nobody".
    assert "users" not in body


def test_authenticator_is_asked_for_the_tas_the_caller_named(
    authenticated_client, monkeypatch, no_ta_service
):
    """The path segment reaches the connector unaltered."""
    import ta_auth_connector

    asked: list[str] = []

    def _record(tas: str):
        asked.append(tas)
        return ta_auth_connector.TasAuthUsersGetResponse(users={"abc12345"})

    monkeypatch.setattr(ta_auth_connector, "get_auth_users", _record)

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    assert asked == [_TAS]
