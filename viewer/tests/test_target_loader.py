"""Seam for tests that need a realistically *loaded* target.

This is the deferred half of the test framework. Today the suite builds the
minimal objects it needs by hand (see ``conftest.py``); the data-heavy model
and API tests will populate the database once instead, by running the real
target loader against a small committed target archive.

To activate this:

1. Commit a small target archive at ``ARCHIVE_PATH`` below
   (the loader's ``.tgz`` bundle format - the historical ``Mpro-v1-zero.tgz``
   used by the old loader test, or a freshly minted minimal one).
2. The ``load_test_target`` fixture and the smoke test below will then run.
3. Build out from there: assert the loader created the expected
   Experiment -> SiteObservation -> Pose graph, then add API/model tests that
   request the fixture for a realistic data set, and port the loader-internals
   tests (``_process_quat_assembly`` etc.) from the old ``test_loader.py``.

Until the archive exists every test in this module is skipped, so the suite
stays green.
"""
# Standard pytest patterns pylint misreads: a test reuses its fixture's name as
# an argument, and django_db_setup is requested for ordering side effect only.
# pylint: disable=redefined-outer-name,unused-argument
from pathlib import Path

import pytest

from viewer.models import Target
from viewer.target_loader import load_target

#: Where a committed test target archive is expected. Does not yet exist.
ARCHIVE_PATH = Path(__file__).parent / "test_data" / "Mpro-v1-zero.tgz"

#: Proposal/visit reference the archive is loaded under (becomes a Project).
TEST_PROPOSAL_REF = "lb-test"

requires_archive = pytest.mark.skipif(
    not ARCHIVE_PATH.is_file(),
    reason=(
        f"No test target archive at {ARCHIVE_PATH}; commit one to enable "
        "loaded-data tests (see module docstring)."
    ),
)


@pytest.fixture(scope="session")
def load_test_target(django_db_setup, django_db_blocker):
    """Load the committed target archive once for the whole test session.

    Session-scoped so the (relatively expensive) load happens a single time;
    tests that need a realistic data graph depend on this fixture. Uses
    ``django_db_blocker`` because session-scoped fixtures sit outside the
    per-test DB transaction that the ``db`` fixture provides.
    """
    if not ARCHIVE_PATH.is_file():
        pytest.skip(f"No test target archive at {ARCHIVE_PATH}")

    with django_db_blocker.unblock():
        load_target(str(ARCHIVE_PATH), proposal_ref=TEST_PROPOSAL_REF)
        yield Target.objects.get(project__title=TEST_PROPOSAL_REF)


@requires_archive
def test_target_loaded(load_test_target):
    """Smoke test: the archive loads and produces a Target. Expand from here."""
    target = load_test_target
    assert target.pk is not None
    assert target.project.title == TEST_PROPOSAL_REF
