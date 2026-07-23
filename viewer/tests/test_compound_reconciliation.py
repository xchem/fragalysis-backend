"""Tests for viewer.compound_reconciliation - the shared compound-matching used
by both the upload validation endpoint and the target loader."""

import pytest

from viewer.compound_reconciliation import (
    MatchStatus,
    inchi_key_for_smiles,
    reconcile_compounds,
)
from viewer.models import Compound, Project

# Two enantiomers - same skeleton, different stereo -> different InChIKey.
R_SMILES = "C[C@H](N)C(=O)O"
S_SMILES = "C[C@@H](N)C(=O)O"
ETHANOL = "CCO"


def _make_project(title="lb00001-1"):
    return Project.objects.create(title=title)


def _make_compound(project, smiles, **fields):
    return Compound.objects.create(
        project=project,
        smiles=smiles,
        inchi=f"InChI={smiles}",
        inchi_key=inchi_key_for_smiles(smiles),
        **fields,
    )


@pytest.mark.django_db
def test_no_existing_match_is_create():
    project = _make_project()
    result = reconcile_compounds(project, [{"smiles": ETHANOL, "compound_code": "X1"}])
    assert [m.status for m in result.matches] == [MatchStatus.CREATE]
    assert not result.needs_curation


@pytest.mark.django_db
def test_project_none_is_all_create():
    result = reconcile_compounds(None, [{"smiles": ETHANOL}])
    assert result.matches[0].status == MatchStatus.CREATE
    assert not result.needs_curation


@pytest.mark.django_db
def test_single_match_identical_is_reuse():
    project = _make_project()
    _make_compound(project, ETHANOL, compound_code="X1", ligand_name="LIG")
    result = reconcile_compounds(
        project, [{"smiles": ETHANOL, "compound_code": "X1", "ligand_name": "LIG"}]
    )
    assert result.matches[0].status == MatchStatus.REUSE
    assert not result.needs_curation


@pytest.mark.django_db
def test_single_match_field_diff_is_conflict():
    project = _make_project()
    _make_compound(project, ETHANOL, compound_code="OLD")
    result = reconcile_compounds(project, [{"smiles": ETHANOL, "compound_code": "NEW"}])
    match = result.matches[0]
    assert match.status == MatchStatus.CONFLICT
    assert match.conflicts["compound_code"] == {"existing": "OLD", "incoming": "NEW"}
    assert result.needs_curation


@pytest.mark.django_db
def test_null_or_blank_is_not_a_conflict():
    project = _make_project()
    # existing has no compound_code; incoming supplies one -> fill, not conflict
    _make_compound(project, ETHANOL, compound_code=None)
    result = reconcile_compounds(project, [{"smiles": ETHANOL, "compound_code": "NEW"}])
    assert result.matches[0].status == MatchStatus.REUSE
    # and the reverse: incoming blank against an existing value is not a conflict
    result2 = reconcile_compounds(project, [{"smiles": ETHANOL, "compound_code": None}])
    assert result2.matches[0].status == MatchStatus.REUSE


@pytest.mark.django_db
def test_multiple_matches_is_ambiguous():
    project = _make_project()
    _make_compound(project, ETHANOL, compound_code="A")
    _make_compound(project, ETHANOL, compound_code="B")
    result = reconcile_compounds(project, [{"smiles": ETHANOL}])
    match = result.matches[0]
    assert match.status == MatchStatus.AMBIGUOUS
    assert len(match.existing) == 2
    assert result.needs_curation


@pytest.mark.django_db
def test_stereo_isomers_do_not_match():
    project = _make_project()
    _make_compound(project, R_SMILES, compound_code="R")
    # incoming is the S enantiomer -> different stereo InChIKey -> no match
    result = reconcile_compounds(project, [{"smiles": S_SMILES}])
    assert result.matches[0].status == MatchStatus.CREATE
    assert inchi_key_for_smiles(R_SMILES) != inchi_key_for_smiles(S_SMILES)


@pytest.mark.django_db
def test_match_is_scoped_to_project():
    project_a = _make_project("lb00001-1")
    project_b = _make_project("lb00002-1")
    _make_compound(project_a, ETHANOL, compound_code="A")
    # same molecule, different project -> no match in project_b
    result = reconcile_compounds(project_b, [{"smiles": ETHANOL}])
    assert result.matches[0].status == MatchStatus.CREATE


@pytest.mark.django_db
def test_unparseable_smiles_is_create():
    project = _make_project()
    result = reconcile_compounds(project, [{"smiles": "not-a-molecule"}])
    assert result.matches[0].status == MatchStatus.CREATE
    assert result.matches[0].inchi_key == ""


@pytest.mark.django_db
def test_curation_payload_shape():
    project = _make_project()
    _make_compound(project, ETHANOL, compound_code="OLD")
    result = reconcile_compounds(project, [{"smiles": ETHANOL, "compound_code": "NEW"}])
    payload = result.curation_payload()
    assert len(payload) == 1
    entry = payload[0]
    assert entry["status"] == "conflict"
    assert entry["conflicts"]["compound_code"] == {"existing": "OLD", "incoming": "NEW"}
    assert entry["incoming"]["compound_code"] == "NEW"
    assert entry["existing"][0]["compound_code"] == "OLD"
