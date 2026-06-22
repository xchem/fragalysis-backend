"""Unit tests for the access-control layer (``api.security.ISPyBSafeQuerySet``).

This is the highest-value area to test: every view that exposes project data
relies on it, and its only external dependency is a single function -
``ta_auth_connector.get_auth_target_access`` - which the
``mock_target_access`` fixture patches.

We exercise the helper methods directly on a bare ``ISPyBSafeQuerySet``
instance (it needs no request for these), rather than through a full view.
"""
# Deliberately reach into _get_q_filter to assert the query it builds.
# pylint: disable=protected-access
from django.db.models import Q

from api.security import ISPyBSafeQuerySet


def test_proposals_from_django_membership(make_user, make_project, settings):
    """With TA_AUTH_SERVICE unset, proposals come from Project.user_id membership."""
    settings.TA_AUTH_SERVICE = ""
    user = make_user("alice")
    make_project("proposal-a", members=[user])
    make_project("proposal-b")  # user is NOT a member

    safe_qs = ISPyBSafeQuerySet()
    proposals = safe_qs.get_proposals_for_user(user, restrict_public_to_membership=True)

    assert proposals == ["proposal-a"]


def test_proposals_from_ta_service(make_user, settings, mock_target_access):
    """With TA_AUTH_SERVICE set, proposals come from the (mocked) TA connector."""
    settings.TA_AUTH_SERVICE = "ta-auth.example"
    user = make_user("bob")
    mock_target_access(["ta-1", "ta-2"])

    safe_qs = ISPyBSafeQuerySet()
    proposals = safe_qs.get_proposals_for_user(user, restrict_public_to_membership=True)

    assert set(proposals) == {"ta-1", "ta-2"}


def test_open_proposals_include_public_projects_and_setting(make_project, settings):
    settings.PUBLIC_TAS_LIST = ["env-public"]
    make_project("open-project", open_to_public=True)
    make_project("closed-project", open_to_public=False)

    safe_qs = ISPyBSafeQuerySet()
    open_proposals = safe_qs.get_open_proposals()

    assert open_proposals == {"open-project", "env-public"}


def test_public_proposals_added_when_not_restricted(make_user, make_project, settings):
    """When restrict_public_to_membership is False, open proposals are folded in."""
    settings.TA_AUTH_SERVICE = ""
    settings.PUBLIC_TAS_LIST = []
    user = make_user("carol")
    make_project("members-only", members=[user])
    make_project("public", open_to_public=True)

    safe_qs = ISPyBSafeQuerySet()
    proposals = safe_qs.get_proposals_for_user(
        user, restrict_public_to_membership=False
    )

    assert set(proposals) == {"members-only", "public"}


def test_user_is_member_of_any_given_proposals(make_user, make_project, settings):
    settings.TA_AUTH_SERVICE = ""
    user = make_user("dave")
    make_project("proposal-x", members=[user])

    safe_qs = ISPyBSafeQuerySet()

    assert safe_qs.user_is_member_of_any_given_proposals(user, ["proposal-x"]) is True
    assert safe_qs.user_is_member_of_any_given_proposals(user, ["proposal-y"]) is False


def test_user_is_member_of_target(make_user, make_project, make_target, settings):
    settings.TA_AUTH_SERVICE = ""
    user = make_user("erin")
    member_target = make_target(make_project("owned", members=[user]))
    other_target = make_target(make_project("not-owned"), title="Other")

    safe_qs = ISPyBSafeQuerySet()

    assert safe_qs.user_is_member_of_target(user, member_target) is True
    assert safe_qs.user_is_member_of_target(user, other_target) is False


def test_get_q_filter_with_permission_string():
    """A filter_permissions string builds a title-in OR open-to-public Q."""
    safe_qs = ISPyBSafeQuerySet()
    safe_qs.filter_permissions = "experiment__experiment_upload__target__project"

    q = safe_qs._get_q_filter(["proposal-a"])

    expected = Q(
        experiment__experiment_upload__target__project__title__in=["proposal-a"]
    ) | Q(experiment__experiment_upload__target__project__open_to_public=True)
    assert q == expected


def test_get_q_filter_without_permission_string():
    """An empty filter_permissions assumes the queryset is over Project itself."""
    safe_qs = ISPyBSafeQuerySet()
    safe_qs.filter_permissions = ""

    q = safe_qs._get_q_filter(["proposal-a"])

    expected = Q(title__in=["proposal-a"]) | Q(open_to_public=True)
    assert q == expected
