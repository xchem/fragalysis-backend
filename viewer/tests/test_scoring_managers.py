"""Tests for the scoring ``by_target`` managers (issue #984, phase 3).

Each scoring model has a ``filter_manager.by_target()`` that prefetches the
SiteObservation->experiment->experiment_upload->target chain and annotates the
target via an ``F`` expression before filtering. These relation paths are easy
to break silently, so this exercises each manager's query end-to-end (it must
compile and run against the schema); an empty target is sufficient to catch a
bad relation path, which raises rather than returning rows.
"""
import pytest

from scoring.models import CmpdChoice, ScoreChoice, SiteObservationChoice, ViewScene


@pytest.mark.django_db
def test_score_choice_by_target_runs(make_project, make_target):
    target = make_target(make_project("proposal-score"))
    assert list(ScoreChoice.filter_manager.by_target(target)) == []


@pytest.mark.django_db
def test_site_observation_choice_by_target_runs(make_project, make_target):
    target = make_target(make_project("proposal-sobs"))
    assert list(SiteObservationChoice.filter_manager.by_target(target)) == []


@pytest.mark.django_db
def test_cmpd_choice_by_target_runs(make_project, make_target):
    target = make_target(make_project("proposal-cmpd"))
    assert list(CmpdChoice.filter_manager.by_target(target)) == []


@pytest.mark.django_db
def test_view_scene_by_target_runs(make_project, make_target):
    target = make_target(make_project("proposal-view"))
    assert list(ViewScene.filter_manager.by_target(target)) == []
