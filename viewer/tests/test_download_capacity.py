"""Tests for the download concurrency cap (ticket #948).

``download_capacity_exceeded`` counts the DownloadLinks rows that are still
being built (a ``task_id`` but no ``file_url``) and compares that against a cap
derived from the celery workers' advertised concurrency. A record that has been
*expired* (``expired_date`` set) is no longer in progress, even if it never got
a ``file_url`` - so it must not be counted (otherwise the cap is reached by dead
records and new downloads are wrongly rejected with HTTP 429).

The single external seam is the celery control inspector, patched here so no
broker/worker is needed.
"""
# pylint: disable=redefined-outer-name,unused-argument
from datetime import datetime, timezone

import pytest

from viewer import download_structures
from viewer.models import DownloadLinks


@pytest.fixture
def now():
    return datetime.now(timezone.utc)


@pytest.fixture
def workers(mocker):
    """Patch celery so inspect().stats() reports a single worker.

    Returns a setter for the worker's pool max-concurrency (default 8).
    """

    def _set(concurrency: int = 8) -> None:
        inspector = mocker.Mock()
        inspector.stats.return_value = {
            "celery@worker1": {"pool": {"max-concurrency": concurrency}}
        }
        mocker.patch.object(
            download_structures.celery_app.control,
            "inspect",
            return_value=inspector,
        )

    _set()
    return _set


@pytest.fixture
def cap_of_one(settings):
    """A percentage that yields a cap of 1 for an 8-wide worker (8 * 1 // 100 -> 1)."""
    settings.MAX_DOWNLOAD_CONCURRENCY_PERCENT = 1


def _make_in_progress(create_date, *, expired_date=None, file_url=None):
    """An in-progress build: a task_id, no file_url (unless given)."""
    return DownloadLinks.objects.create(
        task_id="task-1",
        create_date=create_date,
        file_url=file_url,
        expired_date=expired_date,
    )


@pytest.mark.django_db
def test_disabled_when_percent_is_zero(workers, settings):
    """A percent of 0 disables the guard entirely."""
    settings.MAX_DOWNLOAD_CONCURRENCY_PERCENT = 0
    assert download_structures.download_capacity_exceeded() is False


@pytest.mark.django_db
def test_active_in_progress_record_counts(workers, cap_of_one, now):
    """A genuine in-progress build (no expired_date) is counted."""
    _make_in_progress(now)
    assert download_structures.download_capacity_exceeded() is True


@pytest.mark.django_db
def test_expired_in_progress_records_are_not_counted(workers, cap_of_one, now):
    """The #948 bug: records with an expired_date are no longer in progress and
    must not count towards the cap, even with no file_url."""
    for _ in range(5):
        _make_in_progress(now, expired_date=now)

    assert download_structures.download_capacity_exceeded() is False


@pytest.mark.django_db
def test_completed_record_is_not_counted(workers, cap_of_one, now):
    """A finished build (file_url set) is not in progress."""
    _make_in_progress(now, file_url="Mpro.zip")
    assert download_structures.download_capacity_exceeded() is False


@pytest.mark.django_db
def test_skips_cap_check_when_no_workers_respond(mocker, cap_of_one, now):
    """If no workers respond to stats() the guard is skipped (fail open)."""
    inspector = mocker.Mock()
    inspector.stats.return_value = None
    mocker.patch.object(
        download_structures.celery_app.control, "inspect", return_value=inspector
    )
    _make_in_progress(now)

    assert download_structures.download_capacity_exceeded() is False
