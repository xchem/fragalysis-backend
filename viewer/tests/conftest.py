"""Shared pytest fixtures for the viewer test-suite.

The fixtures here are deliberately *DB-light*: they create only the minimal
objects a unit test needs (a User, a Project, a Target) rather than a full,
loaded data graph. Comprehensive model/API tests that need a realistic data
set will arrive once a real target archive is committed and loaded via the
``load_test_target`` seam (see ``test_target_loader.py``).
"""
# Fixtures that depend on other fixtures legitimately reuse their names as
# arguments, and the `db`/blocker fixtures are requested for their side effect
# only - both are standard pytest patterns that pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument
from typing import Callable, Iterable

import pytest
from django.contrib.auth.models import User
from rest_framework.test import APIClient

import api.ta_auth_connector as ta_auth_connector
from viewer.models import Project, Target


@pytest.fixture
def api_client() -> APIClient:
    """An unauthenticated DRF test client."""
    return APIClient()


@pytest.fixture
def make_user(db) -> Callable[..., User]:
    """Factory returning a freshly created user.

    Usage: ``user = make_user("alice")``.
    """
    counter = {"n": 0}

    def _make_user(username: str | None = None, **kwargs) -> User:
        if username is None:
            counter["n"] += 1
            username = f"user{counter['n']}"
        return User.objects.create_user(username=username, **kwargs)

    return _make_user


@pytest.fixture
def user(make_user) -> User:
    """A single, ready-to-use authenticated-able user."""
    return make_user("tester")


@pytest.fixture
def make_project(db) -> Callable[..., Project]:
    """Factory returning a Project, optionally with members and public flag.

    ``members`` are added to the Project.user_id M2M so the Django-based
    proposal-membership path (ISPyBSafeQuerySet, with TA_AUTH_SERVICE unset)
    treats them as members.
    """

    def _make_project(
        title: str,
        members: Iterable[User] = (),
        open_to_public: bool = False,
    ) -> Project:
        project = Project.objects.create(title=title, open_to_public=open_to_public)
        for member in members:
            project.user_id.add(member)
        return project

    return _make_project


@pytest.fixture
def make_target(db) -> Callable[..., Target]:
    """Factory returning a Target linked to the given project."""

    def _make_target(project: Project, title: str = "Mpro") -> Target:
        return Target.objects.create(title=title, project=project)

    return _make_target


@pytest.fixture
def authenticated_client(api_client, user) -> APIClient:
    """A DRF client authenticated as ``user`` (bypasses real OIDC/session auth)."""
    api_client.force_authenticate(user=user)
    return api_client


@pytest.fixture
def mock_target_access(monkeypatch) -> Callable[[Iterable[str]], None]:
    """Patch the single external access-control seam.

    The security layer (api.security.get_proposals_for_user) calls
    ``ta_auth_connector.get_auth_target_access(username)`` when
    ``settings.TA_AUTH_SERVICE`` is set. This fixture returns a setter that
    makes that call return a chosen set of proposal strings, so tests can drive
    access control through the TA path without any real HTTP service.

    Tests using this must also set ``settings.TA_AUTH_SERVICE`` (e.g. via the
    pytest-django ``settings`` fixture) to enable the TA path.
    """

    def _set(proposals: Iterable[str]) -> None:
        proposal_set = set(proposals)
        monkeypatch.setattr(
            ta_auth_connector,
            "get_auth_target_access",
            lambda username: set(proposal_set),
        )

    return _set
