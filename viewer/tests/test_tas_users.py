"""Tests for the ``/api/tas/`` endpoint (issue #1017).

The endpoint answers "who else is a member of this target access string?", by
way of the authenticator's ``/users/{tas}`` endpoint (TA-Auth 1.5.0). It is
shaped like ``/api/user/``, which answers the mirror-image question, and
carries the same authenticator/ping block. Two things are worth stating up
front, because they shape every test here:

- **It is membership-gated.** A caller may only ask about a TAS they are
  themselves a member of - the answer is a list of people, and knowing who has
  access to a proposal you have no part in is not the caller's business. Public
  (``open_to_public``) proposals do not widen this: the check uses
  ``restrict_public_to_membership=True``, matching the rule that everyone may
  *see* public targets but only members may act on them.
- **"Nobody" is not "we do not know".** The authenticator answers 503, rather
  than an empty set, when it cannot reach ISPyB, and the client preserves that
  distinction. The endpoint must too - reporting an empty membership for a
  service outage would be a silent lie about who has access.

Membership is driven through the Django ``Project.user_id`` M2M (with
``TA_AUTH_SERVICE`` unset), which is the same seam ``test_security_helpers``
uses; the authenticator call itself is mocked via ``mock_tas_users``.
"""
# The `no_ta_service` fixture is requested for its side effect only - the
# standard pytest pattern that pylint misreads (as it does in conftest.py).
# pylint: disable=unused-argument

import pytest

_TAS = "lb12345-1"


@pytest.fixture(name="no_ta_service")
def fixture_no_ta_service(settings):
    """Drive proposal membership from Django rather than the TA service."""
    settings.TA_AUTH_SERVICE = ""


@pytest.mark.django_db
def test_requires_authentication(api_client, mock_tas_users, no_ta_service):
    """An anonymous caller is not told who has access to anything."""
    mock_tas_users(users=["abc12345"])

    response = api_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code in (401, 403)


def test_member_sees_the_users_of_the_tas(
    authenticated_client, user, make_project, mock_tas_users, no_ta_service
):
    """A member gets the membership the authenticator reports."""
    make_project(_TAS, members=[user])
    mock_tas_users(users=["xyz98765", "abc12345"])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    body = response.json()
    assert body["tas"] == _TAS
    # Sorted, so the response is stable for the caller (the source is a set).
    assert body["users"] == ["abc12345", "xyz98765"]
    # The whole contract, so the shape cannot drift away from /api/user/'s.
    assert sorted(body.keys()) == ["authenticator", "ping", "tas", "users"]


def test_response_names_the_authenticator_that_answered(
    authenticated_client,
    user,
    make_project,
    mock_tas_users,
    mock_authenticator,
    no_ta_service,
):
    """The response carries the same authenticator/ping block as /api/user/.

    The two endpoints are a pair - both are answers *from the authenticator* -
    so a caller can tell which service, at which version, produced the list.
    """
    make_project(_TAS, members=[user])
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


def test_non_member_is_refused(
    authenticated_client, make_user, make_project, mock_tas_users, no_ta_service
):
    """A TAS the caller is not a member of is refused, not answered."""
    make_project(_TAS, members=[make_user("someone-else")])
    mock_tas_users(users=["someone-else"])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 403


def test_unknown_tas_is_refused(authenticated_client, mock_tas_users, no_ta_service):
    """A TAS the caller has no relationship with is refused.

    There is no Project at all here - the caller must not be able to probe for
    the existence of proposals by comparing 403 with 404.
    """
    mock_tas_users(users=["abc12345"])

    response = authenticated_client.get("/api/tas/lb99999-9/")

    assert response.status_code == 403


def test_public_proposal_does_not_grant_access(
    authenticated_client, make_project, mock_tas_users, no_ta_service
):
    """Being public makes a proposal *visible*, not its membership readable."""
    make_project(_TAS, open_to_public=True)
    mock_tas_users(users=["abc12345"])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 403


def test_tas_with_no_members_is_an_empty_list_not_an_error(
    authenticated_client, user, make_project, mock_tas_users, no_ta_service
):
    """An empty membership is a legitimate answer, and answered as 200."""
    make_project(_TAS, members=[user])
    mock_tas_users(users=[])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    body = response.json()
    assert body["tas"] == _TAS
    assert body["users"] == []


def test_authenticator_failure_is_reported_not_hidden(
    authenticated_client, user, make_project, mock_tas_users, no_ta_service
):
    """When the authenticator cannot answer, say so - do not report "nobody"."""
    make_project(_TAS, members=[user])
    mock_tas_users(error="Status was 503 (not 200)")

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 503
    body = response.json()
    assert body["tas"] == _TAS
    assert "error" in body
    # No 'users' key at all - an empty list here would read as "nobody".
    assert "users" not in body


def test_membership_is_checked_before_the_authenticator_is_asked(
    authenticated_client, make_project, monkeypatch, no_ta_service
):
    """A refused caller causes no query of the authenticator (or ISPyB).

    '/users/{tas}' is uncached upstream - every call reaches the database - so
    an unauthorised request must be stopped here rather than passed on.
    """
    import ta_auth_connector

    make_project(_TAS)
    calls: list[str] = []

    def _record(tas: str):
        calls.append(tas)
        return ta_auth_connector.TasAuthUsersGetResponse(users=set())

    monkeypatch.setattr(ta_auth_connector, "get_auth_users", _record)

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 403
    assert calls == []


def test_membership_may_come_from_the_ta_service(
    authenticated_client, settings, mock_target_access, mock_tas_users
):
    """Membership is whatever get_proposals_for_user says, TA service included."""
    settings.TA_AUTH_SERVICE = "http://auth.example.svc"
    mock_target_access([_TAS])
    mock_tas_users(users=["abc12345"])

    response = authenticated_client.get(f"/api/tas/{_TAS}/")

    assert response.status_code == 200
    assert response.json()["users"] == ["abc12345"]
