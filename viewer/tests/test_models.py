"""Light model-behaviour tests for the viewer app.

These cover behaviour reachable with a minimal object graph: string
representations and the database-level unique constraints. Heavier behaviour
that needs a fully loaded data graph - e.g.
``SiteObservationQualityStatus.save()`` main-status enforcement
(viewer/models.py) - is deferred to the real-data phase (see
``test_target_loader.py``).
"""
import pytest
from django.db import IntegrityError, transaction

from viewer.models import ComputedSet, Target


def test_project_str(make_project):
    project = make_project("my-proposal")
    assert str(project) == "my-proposal"


def test_target_str(make_project, make_target):
    target = make_target(make_project("proposal"), title="Mpro")
    assert str(target) == "Mpro"


def test_target_unique_per_project(make_project, make_target):
    """(title, project) is unique - a duplicate title in the same project fails."""
    project = make_project("proposal")
    make_target(project, title="Mpro")

    with pytest.raises(IntegrityError):
        with transaction.atomic():
            Target.objects.create(title="Mpro", project=project)


def test_target_same_title_different_project_allowed(make_project, make_target):
    """The same title is fine in a *different* project."""
    target_a = make_target(make_project("proposal-a"), title="Mpro")
    target_b = make_target(make_project("proposal-b"), title="Mpro")

    assert target_a.pk != target_b.pk


def test_computedset_unique_name_per_target(make_project, make_target, user):
    """(name, target) is unique. owner_user is passed explicitly because the
    model default (ANONYMOUS_USER, pk=1) does not exist in a fresh test DB."""
    target = make_target(make_project("proposal"))
    ComputedSet.objects.create(
        name="set-1", target=target, spec_version=1.2, owner_user=user
    )

    with pytest.raises(IntegrityError):
        with transaction.atomic():
            ComputedSet.objects.create(
                name="set-1", target=target, spec_version=1.2, owner_user=user
            )
