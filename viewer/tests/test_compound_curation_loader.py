"""Loader-side compound-curation tests: the reconciliation gate and the
resolution map it builds for process_compound.

The full create/reuse apply runs through the create_objects machinery and is
exercised by the loader integration test; here we cover the gate's decisions
(fail-safe on unresolved, and the reuse/update action map when the user has
supplied a completed curation spreadsheet).
"""

# pylint: disable=protected-access,redefined-outer-name

import logging

import pytest

from viewer.compound_curation import build_curation_xlsx
from viewer.compound_reconciliation import inchi_key_for_smiles, reconcile_compounds
from viewer.models import Compound, Project, SiteObservation
from viewer.target_loader import ProcessedObject, TargetLoader, _row_named_by_pk
from viewer.tests.curation_sheet_helpers import (
    MERGE_ROW,
    edit_group_row,
    existing_row_ids,
)

ETHANOL = "CCO"
CRYSTALS = {
    "Xtal-1": {
        "crystallographic_files": {
            "ligand_cif": {
                "ligands": {"LIG": {"smiles": ETHANOL, "compound_code": "NEW"}}
            }
        }
    }
}


class _Report:
    def __init__(self):
        self.logs = []

    def log(self, level, msg):
        self.logs.append((level, msg))

    def errors(self):
        return [m for lvl, m in self.logs if lvl == logging.ERROR]


def _loader(project, curation_file=None):
    tl = TargetLoader.__new__(TargetLoader)
    tl.project = project
    tl.report = _Report()  # type: ignore[assignment]
    tl.curation_file = curation_file
    tl._compound_resolution = {}
    tl._compound_supersede = {}
    tl._compound_keys = {}
    tl.target = None
    return tl


def _existing(project, code="OLD"):
    return Compound.objects.create(
        project=project,
        smiles=ETHANOL,
        inchi="InChI=CCO",
        inchi_key=inchi_key_for_smiles(ETHANOL),
        compound_code=code,
    )


def test_extract_incoming_compounds_pure():
    assert TargetLoader._extract_incoming_compounds(CRYSTALS) == [
        {"smiles": ETHANOL, "compound_code": "NEW"}
    ]


@pytest.mark.django_db
def test_gate_fails_safe_on_unresolved_conflict():
    project = Project.objects.create(title="lb-1")
    _existing(project)
    tl = _loader(project)

    tl._reconcile_compounds_gate(CRYSTALS)

    assert tl.report.errors()  # ERROR logged -> triggers rollback
    assert tl._compound_resolution == {}  # nothing to reuse without a decision


@pytest.mark.django_db
def test_gate_builds_update_action_from_curation(tmp_path):
    project = Project.objects.create(title="lb-1")
    existing = _existing(project)

    recon = reconcile_compounds(project, [{"smiles": ETHANOL, "compound_code": "NEW"}])
    data = build_curation_xlsx(recon.curation_payload(), target_name="Mpro")

    # user takes the incoming values over the existing row: MERGE on the
    # existing duplicate, with the incoming compound_code picked in the merge row
    existing_id = existing_row_ids(data)[0]
    data = edit_group_row(data, existing_id, action="MERGE")
    data = edit_group_row(data, MERGE_ROW, action="MERGE", compound_code="NEW")

    path = tmp_path / "curation.xlsx"
    path.write_bytes(data)

    tl = _loader(project, curation_file=str(path))
    tl._reconcile_compounds_gate(CRYSTALS)

    assert not tl.report.errors()
    actions = list(tl._compound_resolution.values())
    assert len(actions) == 1
    assert actions[0].op == "update"
    assert actions[0].existing_id == existing.id
    assert actions[0].incoming["compound_code"] == "NEW"


@pytest.mark.django_db
def test_gate_reuse_action_for_exact_match():
    project = Project.objects.create(title="lb-1")
    existing = _existing(project, code="NEW")  # identical -> auto reuse, no conflict
    tl = _loader(project)

    tl._reconcile_compounds_gate(CRYSTALS)

    assert not tl.report.errors()
    actions = list(tl._compound_resolution.values())
    assert len(actions) == 1
    assert actions[0].op == "reuse"
    assert actions[0].existing_id == existing.id


@pytest.mark.django_db
def test_gate_noop_without_incoming():
    project = Project.objects.create(title="lb-1")
    tl = _loader(project)
    tl._reconcile_compounds_gate({})
    assert not tl.report.errors()
    assert tl._compound_resolution == {}


@pytest.mark.django_db
def test_reuse_applies_to_a_compound_no_target_links_to():
    """The bug this guards. The reachability check used
    Compound.filter_manager.by_target(), which finds compounds only through
    their site observations. A compound matched in the project but not yet
    linked to THIS target was therefore invisible, and a KEEP meaning "reuse
    compound N" silently created a duplicate instead - 481 of them in one real
    upload.
    """
    project = Project.objects.create(title="lb-1")
    existing = _existing(project)
    tl = _loader(project)

    # nothing links it to any target, which is precisely the blind spot
    assert not SiteObservation.objects.filter(cmpd=existing).exists()
    assert tl._resolved_compound_applies(existing.pk)


@pytest.mark.django_db
def test_reuse_does_not_reach_outside_the_project():
    """Widening from target to project must not widen further: reconciliation
    matched within one project, so a decision can only ever apply there."""
    project = Project.objects.create(title="lb-1")
    stranger = _existing(Project.objects.create(title="lb-2"))
    tl = _loader(project)

    assert not tl._resolved_compound_applies(stranger.pk)


@pytest.mark.django_db
def test_reuse_action_is_not_downgraded_for_an_unlinked_compound():
    """End of the chain: gate produces a reuse action, and the loader still
    considers it applicable when the compound belongs to no target."""
    project = Project.objects.create(title="lb-1")
    existing = _existing(project, code="NEW")  # exact match -> reuse, no curation
    tl = _loader(project)

    tl._reconcile_compounds_gate(CRYSTALS)

    assert not tl.report.errors()
    actions = list(tl._compound_resolution.values())
    assert [a.op for a in actions] == ["reuse"]
    assert actions[0].existing_id == existing.pk
    assert tl._resolved_compound_applies(actions[0].existing_id)


def _processed(fields):
    return ProcessedObject(
        model_class=Compound, fields=fields, key=("Xtal-1", "LIG"), defaults={}
    )


@pytest.mark.django_db
def test_pk_lookup_finds_a_row_the_target_scoped_manager_cannot_see():
    """Why this matters: create_objects looks instances up through
    by_target(), and on a miss it used to construct the model with whatever
    `fields` held and save it. With fields={"id": N} that is not an INSERT, it
    is an UPDATE of row N that blanks every column absent from `defaults` -- so
    reusing a sibling target's compound wiped the row it was meant to reuse.
    """
    project = Project.objects.create(title="lb-1")
    existing = _existing(project)

    assert not SiteObservation.objects.filter(cmpd=existing).exists()
    assert _row_named_by_pk(_processed({"id": existing.pk})) == existing


@pytest.mark.django_db
def test_pk_lookup_ignores_natural_key_fields():
    """Only an explicit primary key names a row globally; anything else is a
    normal lookup and must still be allowed to fall through to creation."""
    project = Project.objects.create(title="lb-1")
    _existing(project)

    assert _row_named_by_pk(_processed({"compound_code": "OLD"})) is None


@pytest.mark.django_db
def test_pk_lookup_returns_none_when_the_row_really_is_gone():
    project = Project.objects.create(title="lb-1")
    existing = _existing(project)
    missing = existing.pk
    existing.delete()

    assert _row_named_by_pk(_processed({"id": missing})) is None
