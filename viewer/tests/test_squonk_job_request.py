"""Tests for the Squonk2 job-request guards (issue #984, phase 2).

``viewer.squonk_job_request.create_squonk_job`` validates the snapshot's file
transfers and the Squonk2 project before queueing a job. These guards each
raise ``ValueError`` rather than letting a bad request through; they previously
had no tests. We drive the early validation branches with a stub request and
the shared ``mock_squonk_agent`` fixture, so no live Data-Manager API is
contacted. The success path (which queues a real DM job) stays with the
integration suite.
"""
# pylint: disable=redefined-outer-name
from types import SimpleNamespace

import pytest

from viewer import squonk_job_request
from viewer.models import JobFileTransfer, Snapshot
from viewer.squonk2_agent import Squonk2AgentRv


@pytest.fixture
def job_request():
    """Factory for a minimal DRF-like request create_squonk_job understands."""

    def _make(user, snapshot_id):
        return SimpleNamespace(
            session={"oidc_access_token": "token"},
            user=user,
            data={
                "access": "proposal-1",
                "target": 1,
                "snapshot": snapshot_id,
                "session_project": 1,
                "squonk_job_name": "test-job",
                "squonk_job_spec": "{}",
            },
        )

    return _make


@pytest.mark.django_db
def test_missing_file_transfer_raises(job_request, make_user):
    """A snapshot with no JobFileTransfer cannot run a job."""
    request = job_request(make_user("alice"), snapshot_id=999999)

    with pytest.raises(ValueError, match="No JobFileTransfer object"):
        squonk_job_request.create_squonk_job(request)


@pytest.mark.django_db
def test_incomplete_transfer_raises(job_request, make_user):
    """A file transfer that has not reached SUCCESS blocks the job."""
    snapshot = Snapshot.objects.create(title="snap", data="{}")
    JobFileTransfer.objects.create(
        snapshot=snapshot, transfer_status=JobFileTransfer.STARTED
    )
    request = job_request(make_user("bob"), snapshot_id=snapshot.id)

    with pytest.raises(ValueError, match="Job Transfer not complete"):
        squonk_job_request.create_squonk_job(request)


@pytest.mark.django_db
def test_ensure_project_failure_raises(job_request, make_user, mock_squonk_agent):
    """If the agent cannot get/create a Squonk2 Project, the job is refused."""
    snapshot = Snapshot.objects.create(title="snap", data="{}")
    JobFileTransfer.objects.create(
        snapshot=snapshot, transfer_status=JobFileTransfer.SUCCESS
    )
    mock_squonk_agent(ensure_project=Squonk2AgentRv(success=False, msg="no project"))
    request = job_request(make_user("carol"), snapshot_id=snapshot.id)

    with pytest.raises(ValueError, match="failed to get/create a Squonk2 Project"):
        squonk_job_request.create_squonk_job(request)
