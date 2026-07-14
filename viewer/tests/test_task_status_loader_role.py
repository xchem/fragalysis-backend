"""Tests that a ``Loader``-role user can poll ``task_status`` for their uploads.

The upload endpoint (``UploadExperimentUploadView``) grants ``UserRole.LOADER_ROLE``
users a bypass, letting them load data for *any* proposal even without membership
(see ``api.security.check_upload_tas_authorisation``). ``TaskStatusView`` must
grant the same
bypass, otherwise a Loader who uploads to a non-public proposal they are not a
member of receives a ``task_status_url`` they are then forbidden to poll.

The single external seam is the Celery result backend (``viewer.views.AsyncResult``),
patched here so no broker/worker is needed.
"""
# pylint: disable=redefined-outer-name,unused-argument
import uuid
from unittest import mock

import pytest
from django.urls import reverse

from viewer import views
from viewer.models import UserRole


@pytest.fixture
def mock_task_result(monkeypatch):
    """Patch ``viewer.views.AsyncResult`` to a finished task for ``proposal``.

    Returns a stub whose ``info`` carries the ``proposal_ref`` that
    ``TaskStatusView`` reads to resolve the Project and drive access control.
    """

    def _set(proposal: str):
        result = mock.Mock()
        result.state = "SUCCESS"
        result.ready.return_value = True
        result.info = {"proposal_ref": proposal, "description": ["all done"]}
        monkeypatch.setattr(views, "AsyncResult", mock.Mock(return_value=result))
        return result

    return _set


def _status_url() -> str:
    return reverse("viewer:task_status", kwargs={"task_id": uuid.uuid4()})


@pytest.mark.django_db
def test_loader_can_query_task_status_without_membership(
    api_client, make_user, make_project, mock_task_result
):
    """A Loader gets 200 for a non-public proposal they are not a member of."""
    project = make_project("lb00000-1", open_to_public=False)
    mock_task_result(project.title)

    loader = make_user("loader")
    role, _ = UserRole.objects.get_or_create(name=UserRole.LOADER_ROLE)
    role.users.add(loader)
    api_client.force_authenticate(user=loader)

    response = api_client.get(_status_url())

    assert response.status_code == 200, response.content
    assert response.json()["status"] == "SUCCESS"


@pytest.mark.django_db
def test_non_loader_non_member_is_forbidden(
    api_client, make_user, make_project, mock_task_result
):
    """The bypass is scoped to the role: a plain non-member still gets 403."""
    project = make_project("lb00000-1", open_to_public=False)
    mock_task_result(project.title)

    stranger = make_user("stranger")
    api_client.force_authenticate(user=stranger)

    response = api_client.get(_status_url())

    assert response.status_code == 403, response.content
