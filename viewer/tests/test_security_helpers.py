"""Tests for the membership helpers on ``ISPyBSafeQuerySet`` (issue #991).

``user_is_member_of_target`` and ``user_is_member_of_any_given_proposals`` both
sit on top of ``get_proposals_for_user``. With ``TA_AUTH_SERVICE`` unset the
proposal set is derived from the Django ``Project.user_id`` membership M2M, so
``make_project(members=...)`` is enough to drive both the allow and deny paths
without any external access-control service.
"""
from api.security import ISPyBSafeQuerySet


def test_user_is_member_of_target_true(make_user, make_project, make_target, settings):
    settings.TA_AUTH_SERVICE = ""
    user = make_user("member")
    project = make_project("lb12345-1", members=[user])
    target = make_target(project)

    assert ISPyBSafeQuerySet().user_is_member_of_target(user, target) is True


def test_user_is_member_of_target_false(make_user, make_project, make_target, settings):
    settings.TA_AUTH_SERVICE = ""
    member = make_user("member")
    outsider = make_user("outsider")
    project = make_project("lb12345-1", members=[member])
    target = make_target(project)

    assert ISPyBSafeQuerySet().user_is_member_of_target(outsider, target) is False


def test_user_is_member_of_any_given_proposals_true(make_user, make_project, settings):
    settings.TA_AUTH_SERVICE = ""
    user = make_user("member")
    make_project("lb12345-1", members=[user])

    safe_qs = ISPyBSafeQuerySet()
    assert (
        safe_qs.user_is_member_of_any_given_proposals(user, ["lb99999-1", "lb12345-1"])
        is True
    )


def test_user_is_member_of_any_given_proposals_false(make_user, make_project, settings):
    settings.TA_AUTH_SERVICE = ""
    user = make_user("member")
    make_project("lb12345-1", members=[user])

    safe_qs = ISPyBSafeQuerySet()
    assert safe_qs.user_is_member_of_any_given_proposals(user, ["lb99999-1"]) is False
