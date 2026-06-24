"""Tests for the "lost download" housekeeping (ticket #932).

``expire_lost_download_records`` inspects in-progress DownloadLinks records (those
without a ``keep_zip_until``) and expires the ones whose Celery task has crashed or
been lost, recording a human-readable ``expiry_reason``. ``TaskStatusView`` then
surfaces that reason so a user polling a dead download gets a meaningful error.

The single external seam is the Celery result backend, reached via
``viewer.download_structures.AsyncResult`` - patched here so no broker/worker is
needed. The start-grace / max-runtime windows use the settings defaults (4 / 60 min).
"""
# pylint: disable=redefined-outer-name,unused-argument
import uuid
from datetime import datetime, timedelta, timezone

import pytest
from django.urls import reverse

from viewer import download_structures
from viewer.models import DownloadLinks


def _make_link(
    create_date, *, task_id="task-1", keep_zip_until=None, static_link=False
):
    return DownloadLinks.objects.create(
        task_id=task_id,
        create_date=create_date,
        keep_zip_until=keep_zip_until,
        static_link=static_link,
    )


@pytest.fixture
def now():
    return datetime.now(timezone.utc)


@pytest.mark.django_db
def test_within_start_grace_is_left_alone(mock_async_result, now):
    """A young record is given time to start - not even queried."""
    async_result = mock_async_result("PENDING")
    link = _make_link(now - timedelta(minutes=1))

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is None
    assert link.expiry_reason is None
    async_result.assert_not_called()


@pytest.mark.django_db
@pytest.mark.parametrize("state", ["FAILURE", "REVOKED"])
def test_failed_task_is_expired(mock_async_result, now, state):
    mock_async_result(state, info="boom")
    link = _make_link(now - timedelta(minutes=10))

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is not None
    assert link.expiry_reason
    assert "ended unexpectedly" in link.expiry_reason


@pytest.mark.django_db
def test_pending_task_is_expired_as_lost(mock_async_result, now):
    """PENDING past the grace means no worker result - the task is lost."""
    mock_async_result("PENDING")
    link = _make_link(now - timedelta(minutes=10))

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is not None
    assert "lost" in link.expiry_reason.lower()


@pytest.mark.django_db
def test_running_within_max_runtime_is_left_alone(mock_async_result, now):
    """A genuinely long-running download must not be killed prematurely."""
    mock_async_result("PROCESSING")
    link = _make_link(now - timedelta(minutes=10))

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is None
    assert link.expiry_reason is None


@pytest.mark.django_db
def test_running_beyond_max_runtime_is_expired(mock_async_result, now):
    """A 'running' task older than the max runtime is presumed lost (frozen
    result after a worker restart)."""
    mock_async_result("PROCESSING")
    link = _make_link(now - timedelta(minutes=90))

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is not None
    assert "maximum runtime" in link.expiry_reason


@pytest.mark.django_db
def test_no_task_id_is_expired(mock_async_result, now):
    """A record past the grace that never got a task_id was never launched."""
    async_result = mock_async_result("PENDING")
    link = _make_link(now - timedelta(minutes=10), task_id=None)

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is not None
    assert "never launched" in link.expiry_reason
    async_result.assert_not_called()


@pytest.mark.django_db
def test_completed_record_is_not_a_candidate(mock_async_result, now):
    """A record with keep_zip_until set has finished and is ignored here."""
    async_result = mock_async_result("FAILURE")
    link = _make_link(
        now - timedelta(minutes=10), keep_zip_until=now + timedelta(minutes=30)
    )

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is None
    async_result.assert_not_called()


@pytest.mark.django_db
def test_static_link_is_left_alone(mock_async_result, now):
    async_result = mock_async_result("FAILURE")
    link = _make_link(now - timedelta(minutes=10), static_link=True)

    download_structures.expire_lost_download_records()

    link.refresh_from_db()
    assert link.expired_date is None
    async_result.assert_not_called()


@pytest.mark.django_db
def test_soft_erase_ignores_in_progress_records(now):
    """soft_erase no longer time-expires in-progress (keep_zip_until null)
    records - that is now expire_lost_download_records' job."""
    link = _make_link(now - timedelta(hours=2))

    download_structures.soft_erase_out_of_date_download_records()

    link.refresh_from_db()
    assert link.expired_date is None


@pytest.mark.django_db
def test_task_status_view_returns_recorded_reason(api_client, now):
    """A user querying a lost download gets the stored reason as a FAILED status,
    even with no live Celery result."""
    task_id = uuid.uuid4()
    DownloadLinks.objects.create(
        task_id=str(task_id),
        create_date=now,
        expired_date=now,
        expiry_reason="The download task was lost (no worker result)",
    )

    url = reverse("viewer:task_status", kwargs={"task_id": task_id})
    response = api_client.get(url)

    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "FAILED"
    assert body["finished"] is True
    assert body["messages"] == ["The download task was lost (no worker result)"]
