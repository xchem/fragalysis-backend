"""Standalone (no-DB) tests for viewer.compound_curation.

These cover the stable row id and the decision resolver on plain dicts, so they
need neither openpyxl nor Django. The spreadsheet build/parse round trip lives
in test_compound_curation_hierarchical.py.
"""

from viewer.compound_curation import (
    ParsedCuration,
    auto_merge_non_conflicting_compounds,
    compound_row_id,
    needs_curation,
    resolve_curation,
)

# --- row id -----------------------------------------------------------------


def test_row_id_stable_and_distinct():
    a = compound_row_id("AAA", {"smiles": "CCO", "compound_code": "Z1"})
    assert a == compound_row_id("AAA", {"smiles": "CCO", "compound_code": "Z1"})
    assert a != compound_row_id("AAA", {"smiles": "CCO", "compound_code": "Z2"})
    assert a != compound_row_id("BBB", {"smiles": "CCO", "compound_code": "Z1"})


# --- resolution -------------------------------------------------------------


def _match(status, inchi, incoming, existing):
    return {
        "status": status,
        "inchi_key": inchi,
        "incoming": incoming,
        "existing": existing,
        "conflicts": {},
    }


def _conflict_match(inchi="AAA", smiles="CCO", code="Z1", existing_id=5):
    return _match(
        "conflict",
        inchi,
        {"smiles": smiles, "compound_code": code},
        [{"id": existing_id}],
    )


def _decision(
    action, inchi="AAA", smiles="CCO", code="Z1", existing_id=5, incoming=None
):
    """A ParsedCuration as _parse_incoming_molecules_sheet would produce it:
    existing-row decisions are keyed by the *incoming* identity plus the
    existing row's database id, and the incoming row's own action by the bare
    identity. A missing incoming action means CREATE, the sheet's default."""
    # Heterogeneous keys by design: 4-tuple per existing row, 3-tuple for the
    # incoming row's own action.
    field_actions: dict = {(inchi, smiles, code, existing_id): action}
    if incoming is not None:
        field_actions[(inchi, smiles, code)] = incoming
    return ParsedCuration(payload_hash="", field_actions=field_actions)


def test_resolve_create_and_reuse_need_no_decision():
    matches = [
        _match("create", "K1", {"smiles": "C"}, []),
        _match("reuse", "K2", {"smiles": "CC"}, [{"id": 7}]),
    ]
    plan = resolve_curation(matches, None)
    assert plan.ok
    assert plan.actions[0].op == "create"
    assert plan.actions[1].op == "reuse" and plan.actions[1].existing_id == 7


def test_resolve_conflict_keep_plus_create_adds_the_incoming_compound():
    """The sheet's defaults, and what they say on the tin: the existing compound
    is kept as it is and the incoming one is added alongside."""
    matches = [_conflict_match()]
    plan = resolve_curation(matches, _decision("KEEP"))
    assert plan.ok
    assert plan.actions[0].op == "create" and plan.actions[0].existing_id is None


def test_resolve_conflict_keep_becomes_reuse():
    """Dropping CREATE from the incoming row folds it into the kept row."""
    matches = [_conflict_match()]
    plan = resolve_curation(matches, _decision("KEEP", incoming="KEEP"))
    assert plan.ok
    assert plan.actions[0].op == "reuse" and plan.actions[0].existing_id == 5


def test_resolve_conflict_delete_becomes_create():
    matches = [_conflict_match()]
    plan = resolve_curation(matches, _decision("DELETE"))
    assert plan.ok
    assert plan.actions[0].op == "create"
    assert plan.actions[0].existing_id is None


def test_resolve_conflict_merge_becomes_update_with_chosen_values():
    matches = [_conflict_match()]
    parsed = _decision("MERGE")
    parsed.field_actions[("MERGE", ("AAA", "CCO", "Z1"))] = {"ligand_name": "PICKED"}

    plan = resolve_curation(matches, parsed)

    assert plan.ok
    assert plan.actions[0].op == "update" and plan.actions[0].existing_id == 5
    assert plan.actions[0].incoming["ligand_name"] == "PICKED"


# --- auto-merge ---------------------------------------------------------------


def _null_vs_value_conflict():
    """A conflict auto-merge can settle alone: the existing row has no
    ligand_name, the incoming one does. Nothing for a user to decide."""
    return [
        {
            "status": "conflict",
            "inchi_key": "AAA",
            "incoming": {"smiles": "CCO", "compound_code": "NEW", "ligand_name": "LIG"},
            "existing": [
                {
                    "id": 7,
                    "smiles": "CCO",
                    "compound_code": "NEW",
                    "ligand_name": None,
                }
            ],
            "conflicts": {"ligand_name": {"existing": None, "incoming": "LIG"}},
        }
    ]


def test_auto_merge_keeps_the_entry_rather_than_dropping_it():
    payload = _null_vs_value_conflict()
    out, stats = auto_merge_non_conflicting_compounds(payload)

    assert stats["auto_merged"] == 1
    # same length and order as the input -- callers index the two together
    assert len(out) == len(payload)
    assert out[0]["status"] == "auto_merged"
    assert out[0]["incoming"]["ligand_name"] == "LIG"


def test_auto_merged_entries_are_not_offered_for_curation():
    out, _ = auto_merge_non_conflicting_compounds(_null_vs_value_conflict())
    assert needs_curation(out) == []


def test_auto_merged_compound_reuses_the_existing_row():
    """The regression that mattered: a dropped entry produced no action, so the
    loader created a duplicate compound instead of updating the match."""
    out, _ = auto_merge_non_conflicting_compounds(_null_vs_value_conflict())
    plan = resolve_curation(out, None)

    assert plan.ok
    (action,) = plan.actions
    assert action.op == "update"
    assert action.existing_id == 7
    assert action.incoming["ligand_name"] == "LIG"


def test_auto_merge_leaves_genuine_conflicts_for_the_user():
    payload = _null_vs_value_conflict()
    payload[0]["existing"][0]["ligand_name"] = "OTHER"
    payload[0]["conflicts"]["ligand_name"]["existing"] = "OTHER"

    out, stats = auto_merge_non_conflicting_compounds(payload)

    assert stats["auto_merged"] == 0
    assert stats["genuine_conflicts"] == 1
    assert len(needs_curation(out)) == 1
    # and with no decision supplied it must fail safe, not guess
    assert not resolve_curation(out, None).ok


def test_resolve_unresolved_conflict_without_decision_fails_safe():
    matches = [_conflict_match()]
    plan = resolve_curation(matches, None)
    assert not plan.ok
    assert plan.unresolved[0].inchi_key == "AAA"
