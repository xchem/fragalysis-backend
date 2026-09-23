"""Throwaway smoke exercise for viewer.pose_curation - NOT the real test suite.

Plants one of each defect and drives every finder and planner, so the ORM paths
are proven against a real database rather than merely linted.
"""
# pytest reuses a fixture's name as the test argument by design.
# pylint: disable=redefined-outer-name,unused-argument
from datetime import datetime, timezone
from typing import Any, Callable

import pytest

from viewer.models import (
    CanonSite,
    CanonSiteConf,
    Compound,
    Experiment,
    ExperimentUpload,
    Pose,
    SiteObservation,
    Xtalform,
    XtalformSite,
)
from viewer.pose_curation import (
    NO_SUCCESSOR_DEAD,
    SUCCESSOR_SAME_POSE,
    diagnose_superseded_mains,
    find_anchor_violations,
    find_dead_poses,
    find_display_name_mismatches,
    find_duplicate_live_identities,
    find_main_not_in_pose,
    find_orphan_observations,
    find_poses_with_superseded_main,
    find_poses_without_main,
    find_split_groups,
    find_successors,
    find_superseded_members,
    find_version_inversions,
    find_wrongly_superseded,
    plan_delete_dead_poses,
    plan_fix_superseded_mains,
    plan_sever_superseded_members,
    plan_unsupersede,
    survey,
)


@pytest.fixture
def graph(db, user, make_project, make_target):
    project = make_project("proposal", members=[user])
    target = make_target(project, title="Peel")

    xtalform = Xtalform.objects.create(name="xf1")
    canon_site = CanonSite.objects.create(
        name="cs1", residues=[], version=1, canon_site_num=1
    )
    conf_site = CanonSiteConf.objects.create(
        name="csc1", canon_site=canon_site, residues=[], version=1
    )
    xtalform_site = XtalformSite.objects.create(
        xtalform_site_id="xs1",
        xtalform=xtalform,
        canon_site=canon_site,
        lig_chain="A",
        residues=[],
        version=1,
    )
    compound = Compound.objects.create(
        smiles="CCO", inchi="InChI=1S/C2H6O", project=project
    )
    upload = ExperimentUpload.objects.create(
        project=project,
        target=target,
        committer=user,
        commit_datetime=datetime(2026, 1, 1, tzinfo=timezone.utc),
        upload_data_dir="upload_1",
        upload_version=1,
        data_version_major=3,
        data_version_minor=1,
    )
    experiment = Experiment.objects.create(
        experiment_upload=upload, code="Peel-x0001", xtalform=xtalform
    )

    def obs(code, version, superseded, pose=None, **kw):
        return SiteObservation.objects.create(
            code=code,
            longcode=f"{code}_long",
            experiment=kw.pop("experiment", experiment),
            cmpd=kw.pop("cmpd", compound),
            xtalform_site=xtalform_site,
            canon_site_conf=kw.pop("canon_site_conf", conf_site),
            pose=pose,
            seq_id=kw.pop("seq_id", 501),
            chain_id="A",
            altloc="0",
            smiles="CCO",
            version=version,
            superseded=superseded,
            **kw,
        )

    # the target-scoping path: canon site -> ref conf -> ref observation -> target
    # distinct seq_id: the anchor must NOT share v1/v2's supersede identity
    anchor = obs("Peel-anchor", 1, False, seq_id=999)
    conf_site.ref_site_observation = anchor
    conf_site.save()
    canon_site.ref_conf_site = conf_site
    canon_site.save()

    # DEFECT: pose whose main is superseded, successor is already a member
    pose = Pose.objects.create(
        canon_site=canon_site, compound=compound, display_name="Peel-x0001_v1"
    )
    v1 = obs("Peel-x0001_v1", 1, True, pose=pose)
    v2 = obs("Peel-x0001_v2", 2, False, pose=pose)
    pose.main_site_observation = v1
    pose.save()

    return {
        "target": target,
        "pose": pose,
        "v1": v1,
        "v2": v2,
        "anchor": anchor,
        "compound": compound,
        "canon_site": canon_site,
    }


def test_survey_runs_both_scoped_and_unscoped(graph):
    for scope in (None, "Peel", graph["target"]):
        s = survey(scope)
        assert s.counts["poses with superseded main"] == 1
        assert s.counts["superseded members (rows)"] == 1
    print("\n" + str(survey("Peel")))


def test_every_finder_compiles(graph):
    """The point: prove each ORM path resolves against a real database."""
    finders: list[Callable[..., Any]] = [
        find_poses_with_superseded_main,
        find_superseded_members,
        find_main_not_in_pose,
        find_poses_without_main,
        find_dead_poses,
        find_orphan_observations,
        find_split_groups,
        find_anchor_violations,
        find_duplicate_live_identities,
        find_version_inversions,
        find_wrongly_superseded,
    ]
    for finder in finders:
        assert list(finder("Peel")) is not None, finder.__name__
        assert list(finder(None)) is not None, finder.__name__
    assert find_display_name_mismatches("Peel") is not None


def test_diagnose_and_fix_superseded_main(graph):
    diags = diagnose_superseded_mains("Peel")
    assert len(diags) == 1
    assert diags[0].kind == SUCCESSOR_SAME_POSE
    assert diags[0].successor is not None
    assert diags[0].successor.pk == graph["v2"].pk
    print("\n" + diags[0].explain())

    plan = plan_fix_superseded_mains("Peel")
    print("\n" + str(plan))
    assert len(plan.actions) == 1
    plan.apply()

    graph["pose"].refresh_from_db()
    assert graph["pose"].main_site_observation_id == graph["v2"].pk


def test_sever_skips_the_main_then_severs_it(graph):
    # v1 is still the main, so severing must refuse
    plan = plan_sever_superseded_members("Peel")
    assert plan.actions == []
    assert len(plan.skipped) == 1

    plan_fix_superseded_mains("Peel").apply()

    plan = plan_sever_superseded_members("Peel")
    assert len(plan.actions) == 1
    plan.apply()
    graph["v1"].refresh_from_db()
    assert graph["v1"].pose_id is None


def test_dead_pose_is_diagnosed_not_fixed(graph):
    # make the pose dead: no live members, main superseded, no successor
    graph["v2"].delete()
    diags = diagnose_superseded_mains("Peel")
    assert diags[0].kind == NO_SUCCESSOR_DEAD
    plan = plan_fix_superseded_mains("Peel")
    assert plan.actions == []
    assert plan.skipped

    destructive = plan_delete_dead_poses("Peel")
    assert destructive.destructive
    with pytest.raises(RuntimeError):
        destructive.apply()


def test_null_identity_is_refused(graph):
    orphan = SiteObservation.objects.create(
        code="orphan", longcode="orphan_long", smiles="C", superseded=True
    )
    assert find_successors(orphan).count() == 0
    plan = plan_unsupersede([orphan.pk])
    assert len(plan.actions) == 1


def test_unsupersede_refuses_when_a_live_twin_exists(graph):
    plan = plan_unsupersede([graph["v1"].pk])
    assert plan.actions == []
    assert "live row" in plan.skipped[0][1]
