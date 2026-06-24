"""Tests for the ``RESTRICTED_TAS_USERS`` proposal-membership simulation.

In a real deployment a user's accessible proposals ("Target Access Strings",
e.g. ``lb12345-1``) come from Keycloak plus the target-access authoriser. For
non-production deployments the ``RESTRICTED_TAS_USERS`` environment variable
provides a debug seam that arbitrarily grants a named user a proposal, letting
tests stand in for that whole chain (see ``api.security`` and the
``set_restricted_tas_users`` fixture in ``conftest``).

These tests verify the seam itself, how it feeds the security layer, and that
it surfaces end-to-end through a real access-controlled endpoint - all without
needing a single ``Target`` in the database.
"""

from api.security import ISPyBSafeQuerySet, get_restricted_tas_user_proposal


def test_restricted_tas_user_matches(make_user, set_restricted_tas_users):
    """A user named in RESTRICTED_TAS_USERS is granted the paired proposal."""
    user = make_user("user-a")
    set_restricted_tas_users("user-a:lb12345-1")

    assert get_restricted_tas_user_proposal(user) == {"lb12345-1"}


def test_restricted_tas_user_no_match(make_user, set_restricted_tas_users):
    """A user absent from RESTRICTED_TAS_USERS is granted nothing."""
    user = make_user("user-b")
    set_restricted_tas_users("user-a:lb12345-1")

    assert get_restricted_tas_user_proposal(user) == set()


def test_restricted_tas_user_multiple_entries(make_user, set_restricted_tas_users):
    """Only the entries for the requested user are collected, and they accrue."""
    user = make_user("user-a")
    set_restricted_tas_users("user-a:lb12345-1,user-b:lb99999-1,user-a:lb55555-1")

    assert get_restricted_tas_user_proposal(user) == {"lb12345-1", "lb55555-1"}


def test_restricted_tas_user_empty_setting(make_user, set_restricted_tas_users):
    """An unset variable grants nothing (and must not raise)."""
    user = make_user("user-a")
    set_restricted_tas_users("")

    assert get_restricted_tas_user_proposal(user) == set()


def test_restricted_tas_users_ignored_in_production(
    make_user, settings, set_restricted_tas_users
):
    """The debug seam is inert in production, regardless of the variable."""
    settings.DEPLOYMENT_MODE = "PRODUCTION"
    user = make_user("user-a")
    set_restricted_tas_users("user-a:lb12345-1")

    assert get_restricted_tas_user_proposal(user) == set()


def test_restricted_tas_feeds_proposals_for_user(
    make_user, settings, set_restricted_tas_users
):
    """The granted proposal reaches get_proposals_for_user (the security layer)."""
    settings.TA_AUTH_SERVICE = ""
    user = make_user("user-a")
    set_restricted_tas_users("user-a:lb12345-1")

    safe_qs = ISPyBSafeQuerySet()
    proposals = safe_qs.get_proposals_for_user(user, restrict_public_to_membership=True)

    assert proposals == ["lb12345-1"]


def test_restricted_tas_user_sees_only_granted_proposal(
    api_client, make_user, make_project, settings, set_restricted_tas_users
):
    """End-to-end: the granted user sees their proposal's Project and no other.

    ``/api/projects/`` is access-controlled by ISPyBSafeQuerySet, so this proves
    the simulated membership flows all the way through a real request - with no
    Target involved.
    """
    settings.TA_AUTH_SERVICE = ""
    make_project("lb12345-1")  # the granted proposal (no Target attached)
    make_project("lb99999-1")  # an unrelated proposal
    set_restricted_tas_users("user-a:lb12345-1")

    user_a = make_user("user-a")
    api_client.force_authenticate(user=user_a)
    response = api_client.get("/api/projects/")

    assert response.status_code == 200
    # The Project model's 'title' is exposed as 'target_access_string'.
    proposals = {row["target_access_string"] for row in response.data["results"]}
    assert proposals == {"lb12345-1"}


def test_user_without_grant_sees_no_proposals(
    api_client, make_user, make_project, settings, set_restricted_tas_users
):
    """A different user, not named in the variable, sees no (non-public) proposal."""
    settings.TA_AUTH_SERVICE = ""
    make_project("lb12345-1")
    set_restricted_tas_users("user-a:lb12345-1")

    user_b = make_user("user-b")
    api_client.force_authenticate(user=user_b)
    response = api_client.get("/api/projects/")

    assert response.status_code == 200
    assert response.data["results"] == []
