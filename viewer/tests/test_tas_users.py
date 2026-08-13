"""Tests for the ``/api/tas/`` endpoint (issue #1017).

The endpoint answers "who is a member of this target access string?", by way of
the authenticator's ``/users/{tas}`` endpoint (TA-Auth 1.5.0). It is shaped like
``/api/user/``, which answers the mirror-image question, and carries the same
authenticator/ping block. Four things are worth stating up front, because they
shape every test here:

- **The TAS arrives as a ``?tas=`` query parameter**, which makes this a DRF
  *list* route. That is what gets the endpoint listed in the browsable API root
  - ``APIRootView`` builds the index by reversing each viewset's ``-list`` route
  and silently skips the ones that fail, so a detail-only route
  (``/api/tas/<tas>/``) is undiscoverable there.
- **Authentication is the only gate.** Any logged-in user may ask about any TAS,
  whether or not they are a member of it. There is no membership check - the
  tests below assert that a caller with no relationship to the TAS is answered,
  not refused. Only an anonymous caller is turned away.
- **A malformed TAS is the caller's error, not an outage.** The string is
  checked against ``settings.TAS_REGEX`` (via ``api.utils.validate_tas``, the
  same gate the upload serializers use) before the authenticator is contacted,
  so nonsense comes back as 400 rather than being reported as a 503.
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

    response = api_client.get(f"/api/tas/?tas={_TAS}")

    assert response.status_code in (401, 403)


def test_authenticated_user_sees_the_users_of_the_tas(
    authenticated_client, mock_tas_users, no_ta_service
):
    """A logged-in caller gets the membership the authenticator reports."""
    mock_tas_users(users=["xyz98765", "abc12345"])

    response = authenticated_client.get(f"/api/tas/?tas={_TAS}")

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

    response = authenticated_client.get(f"/api/tas/?tas={_TAS}")

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

    response = authenticated_client.get("/api/tas/?tas=lb99999-9")

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

    response = authenticated_client.get(f"/api/tas/?tas={_TAS}")

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

    response = authenticated_client.get(f"/api/tas/?tas={_TAS}")

    assert response.status_code == 200
    body = response.json()
    assert body["tas"] == _TAS
    assert body["users"] == []


def test_authenticator_failure_is_reported_not_hidden(
    authenticated_client, mock_tas_users, no_ta_service
):
    """When the authenticator cannot answer, say so - do not report "nobody"."""
    mock_tas_users(error="Status was 503 (not 200)")

    response = authenticated_client.get(f"/api/tas/?tas={_TAS}")

    assert response.status_code == 503
    body = response.json()
    assert body["tas"] == _TAS
    assert "error" in body
    # No 'users' key at all - an empty list here would read as "nobody".
    assert "users" not in body


def test_endpoint_is_listed_in_the_api_root(authenticated_client):
    """The endpoint is discoverable from /api/.

    This is the whole reason the TAS is a query parameter rather than a path
    segment: ``APIRootView`` indexes viewsets by reversing '<basename>-list',
    and quietly drops any that has no list route.
    """
    response = authenticated_client.get("/api/")

    assert response.status_code == 200
    assert "tas" in response.json()


def test_missing_tas_parameter_is_a_400(
    authenticated_client, mock_tas_users, no_ta_service
):
    """/api/tas/ with no ?tas= says what it wants, rather than 500-ing.

    The API root links here without a parameter, so this is the first thing a
    caller browsing the API will see - it needs to be useful.
    """
    mock_tas_users(users=["abc12345"])

    response = authenticated_client.get("/api/tas/")

    assert response.status_code == 400
    body = response.json()
    assert "error" in body
    assert "tas" in body["error"]
    assert "users" not in body


@pytest.mark.parametrize(
    "bad_tas",
    [
        "xx12345-1",  # wrong prefix
        "lb1234-1",  # too few digits
        "lb12345",  # missing visit
        "not-a-tas",
    ],
)
def test_malformed_tas_is_a_400(
    authenticated_client, mock_tas_users, no_ta_service, bad_tas
):
    """A string that is not a TAS is the caller's mistake - say so."""
    mock_tas_users(users=["abc12345"])

    response = authenticated_client.get(f"/api/tas/?tas={bad_tas}")

    assert response.status_code == 400
    body = response.json()
    assert body["tas"] == bad_tas
    assert "error" in body
    assert "users" not in body


def test_malformed_tas_is_not_passed_to_the_authenticator(
    authenticated_client, monkeypatch, no_ta_service
):
    """Validation happens first, so a bad TAS costs no ISPyB query.

    Without this the authenticator would answer 400, which the client reports
    as a bare error string - indistinguishable from an outage, and returned to
    the caller as a misleading 503.
    """
    import ta_auth_connector

    asked: list[str] = []

    def _record(tas: str):
        asked.append(tas)
        return ta_auth_connector.TasAuthUsersGetResponse(users=set())

    monkeypatch.setattr(ta_auth_connector, "get_auth_users", _record)

    response = authenticated_client.get("/api/tas/?tas=not-a-tas")

    assert response.status_code == 400
    assert asked == []


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

    response = authenticated_client.get(f"/api/tas/?tas={_TAS}")

    assert response.status_code == 200
    assert asked == [_TAS]
