"""Tests for ``viewer.permissions.IsObjectProposalMember`` (issue #991).

The permission grants object-level write access only to users who are members of
one of the object's proposals. It resolves the object's proposals via the view's
``filter_permissions`` attribute (a relation path like ``project``). These tests
drive it with a real ``Target`` object and a ``RequestFactory`` request, using
``make_project(members=...)`` to control membership - no external service.
"""
# pylint: disable=redefined-outer-name
import pytest
from django.contrib.auth.models import AnonymousUser
from django.test import RequestFactory
from rest_framework.exceptions import PermissionDenied

from viewer.permissions import IsObjectProposalMember


class _View:
    """A stand-in DRF view exposing only what the permission reads."""

    def __init__(self, filter_permissions="project"):
        self.filter_permissions = filter_permissions


@pytest.fixture
def request_factory():
    return RequestFactory()


def test_member_is_granted(
    request_factory, make_user, make_project, make_target, settings
):
    settings.TA_AUTH_SERVICE = ""
    user = make_user("member")
    target = make_target(make_project("lb12345-1", members=[user]))

    request = request_factory.get("/")
    request.user = user

    permission = IsObjectProposalMember()
    assert permission.has_object_permission(request, _View(), target) is True


def test_non_member_is_denied(
    request_factory, make_user, make_project, make_target, settings
):
    settings.TA_AUTH_SERVICE = ""
    member = make_user("member")
    outsider = make_user("outsider")
    target = make_target(make_project("lb12345-1", members=[member]))

    request = request_factory.get("/")
    request.user = outsider

    permission = IsObjectProposalMember()
    with pytest.raises(PermissionDenied):
        permission.has_object_permission(request, _View(), target)


def test_unauthenticated_user_is_denied(request_factory, make_project, make_target):
    target = make_target(make_project("lb12345-1"))

    request = request_factory.get("/")
    request.user = AnonymousUser()

    permission = IsObjectProposalMember()
    assert permission.has_object_permission(request, _View(), target) is False


def test_view_without_filter_permissions_raises(
    request_factory, make_user, make_project, make_target
):
    user = make_user("member")
    target = make_target(make_project("lb12345-1", members=[user]))

    request = request_factory.get("/")
    request.user = user

    class _BadView:
        pass

    permission = IsObjectProposalMember()
    with pytest.raises(AttributeError):
        permission.has_object_permission(request, _BadView(), target)
