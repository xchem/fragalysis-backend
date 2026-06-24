"""Regression test for ticket #939.

SoakDB-sourced ``FloatField``s on ``Experiment`` can hold ``NaN`` (Postgres
stores it happily). DRF's JSON renderer runs in strict mode (``allow_nan=False``)
and raises ``ValueError: Out of range float values are not JSON compliant: nan``
when rendering, which surfaced as a 500 on ``/api/experiments/``.

``ExperimentReadSerializer`` now maps non-finite floats to ``null`` so the
endpoint renders valid JSON instead of erroring.
"""
# Fixtures legitimately reuse their names as arguments, and the `db` fixture is
# requested for its side effect only - both are standard pytest patterns that
# pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument
from datetime import datetime
from datetime import timezone as dt_timezone

import pytest

from viewer.models import Experiment, ExperimentUpload


@pytest.fixture
def make_experiment(db):
    """Factory creating an Experiment (and its required ExperimentUpload)."""

    def _make_experiment(project, target, committer, **fields):
        upload = ExperimentUpload.objects.create(
            project=project,
            target=target,
            file="experiment-upload/dummy.tgz",
            commit_datetime=datetime(2026, 1, 1, tzinfo=dt_timezone.utc),
            committer=committer,
        )
        return Experiment.objects.create(experiment_upload=upload, **fields)

    return _make_experiment


def test_experiments_list_renders_nan_float_as_null(
    authenticated_client, user, make_project, make_target, make_experiment
):
    """A NaN float field renders as JSON null, not a 500."""
    project = make_project("members-only", members=[user])
    target = make_target(project)
    make_experiment(
        project,
        target,
        committer=user,
        code="x0001",
        dimple_rfree=float("nan"),
    )

    response = authenticated_client.get("/api/experiments/")

    assert response.status_code == 200
    # render the payload to confirm strict JSON encoding succeeds end-to-end
    response.render()
    (row,) = response.data["results"]
    assert row["dimple_rfree"] is None
