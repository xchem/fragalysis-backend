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
from contextlib import contextmanager
from typing import Callable, Iterable
from unittest import mock

import pytest
import ta_auth_connector
from django.contrib.auth.models import User
from rest_framework.test import APIClient

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
def set_restricted_tas_users(settings) -> Callable[[str], None]:
    """Simulate proposal membership via the RESTRICTED_TAS_USERS debug seam.

    ``api.security.get_restricted_tas_user_proposal`` reads two settings that
    the application derives from a single ``RESTRICTED_TAS_USERS`` environment
    variable (see ``fragalysis.settings``): the raw string and its comma-split
    list form. This fixture returns a setter that takes the same
    ``"user:tas,user:tas"`` value the env var uses and applies it to both, so a
    test can grant a user access to a proposal without any external service.
    """

    def _set(value: str) -> None:
        settings.RESTRICTED_TAS_USERS = value
        settings.RESTRICTED_TAS_USERS_LIST = value.split(",") if value else []

    return _set


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


class _FakeNeo4jSession:
    """A neo4j session whose ``run`` raises a chosen exception, or no-ops.

    With ``exc=None`` the session behaves as a reachable graph (``run``
    returns ``None``); pass an exception to simulate an unreachable graph.
    """

    def __init__(self, exc: Exception | None):
        self._exc = exc

    def run(self, *_args, **_kwargs):
        if self._exc is not None:
            raise self._exc
        return None


class _FakeNeo4jDriver:
    """A neo4j driver handing out :class:`_FakeNeo4jSession` instances."""

    def __init__(self, exc: Exception | None):
        self._exc = exc

    @contextmanager
    def session(self):
        yield _FakeNeo4jSession(self._exc)


@pytest.fixture
def mock_neo4j(monkeypatch) -> Callable[..., None]:
    """Patch the neo4j driver seam used by the service-status probes.

    ``service_status.services.get_driver`` is the single point where a real
    neo4j connection is opened. This fixture returns a setter that replaces it
    with a fake driver, so probe tests never touch a real graph. Call with no
    argument for a reachable graph, or pass an exception (e.g.
    ``neo4j.exceptions.ServiceUnavailable(...)``) to simulate an outage.
    """

    def _set(exc: Exception | None = None) -> None:
        from service_status import services

        monkeypatch.setattr(
            services, "get_driver", lambda **_kwargs: _FakeNeo4jDriver(exc)
        )

    return _set


@pytest.fixture
def mock_async_result(monkeypatch) -> Callable[..., mock.Mock]:
    """Patch ``viewer.download_structures.AsyncResult`` to a fixed state.

    The download housekeeping/capacity logic inspects a Celery
    ``AsyncResult``'s ``state`` and ``info`` to decide whether a download task
    is still running, lost, or failed. This fixture returns a setter that makes
    every ``AsyncResult(...)`` return a stub with the given ``state``/``info``,
    so tests need no broker or worker.

    The *constructor* mock is returned (mirroring ``mocker.patch.object``), so a
    test can assert whether the result backend was consulted at all (e.g.
    ``async_result.assert_not_called()``); the stub result is reachable via its
    ``.return_value``.
    """

    def _set(state: str, info=None) -> mock.Mock:
        from viewer import download_structures

        result = mock.Mock()
        result.state = state
        result.info = info
        constructor = mock.Mock(return_value=result)
        monkeypatch.setattr(download_structures, "AsyncResult", constructor)
        return constructor

    return _set


@pytest.fixture
def mock_squonk_agent(monkeypatch) -> Callable[..., mock.Mock]:
    """Patch the Squonk2 agent singleton used by ``squonk_job_request``.

    ``viewer.squonk_job_request`` binds a module-level ``Squonk2Agent``
    singleton (``_SQ2A``) at import time. This fixture returns a setter that
    swaps it for a ``Mock`` whose ``ping``/``ensure_project`` (etc.) return a
    ``Squonk2AgentRv``-shaped value, so the validation/error branches can be
    exercised without a live Data-Manager/Account-Server API. The stub agent
    mock is returned so a test can tailor individual method return values.
    """

    def _set(**method_return_values) -> mock.Mock:
        from viewer import squonk_job_request
        from viewer.squonk2_agent import SuccessRv

        agent = mock.Mock()
        # Default every called method to "success" unless overridden.
        agent.ping.return_value = SuccessRv
        agent.ensure_project.return_value = SuccessRv
        for name, value in method_return_values.items():
            getattr(agent, name).return_value = value
        monkeypatch.setattr(squonk_job_request, "_SQ2A", agent)
        return agent

    return _set
