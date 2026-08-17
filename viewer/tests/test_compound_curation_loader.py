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
from viewer.models import Compound, Project
from viewer.target_loader import TargetLoader
from viewer.tests.curation_sheet_helpers import edit_group_row, existing_row_ids

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
    data = edit_group_row(data, "Merge", action="MERGE", compound_code="NEW")

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
