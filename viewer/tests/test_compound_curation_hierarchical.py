"""Standalone (no-DB) tests for the hierarchical compound-group curation sheet.

These cover the build -> edit -> parse round trip of the "Incoming compounds"
sheet, and in particular that its column layout is derived from CONTENT_FIELDS
on both sides rather than hardcoded. They need only openpyxl.
"""

import io

import pytest
from openpyxl import load_workbook

from viewer import compound_curation as cc
from viewer.compound_curation import (
    CurationFormatError,
    build_curation_xlsx,
    curation_headers,
    parse_curation_xlsx,
    resolve_curation,
)
from viewer.tests.curation_sheet_helpers import SHEET
from viewer.tests.curation_sheet_helpers import edit_group_row as _edit
from viewer.tests.curation_sheet_helpers import group_headers as _headers_in

INCOMING = {
    "smiles": "CCO",
    "compound_code": "NEW-1",
    "ligand_name": "LIG-NEW",
    "modeled_smiles_canon": "CCO",
}
EXISTING = {
    "id": 42,
    "inchi_key": "AAA",
    "smiles": "CCO",
    "compound_code": "OLD-1",
    "ligand_name": "LIG-OLD",
    "modeled_smiles_canon": "CCO",
}
KEY = ("AAA", "CCO", "NEW-1")


def _conflict_payload():
    return [
        {
            "status": "conflict",
            "inchi_key": "AAA",
            "incoming": dict(INCOMING),
            "conflicts": {
                "compound_code": {"existing": "OLD-1", "incoming": "NEW-1"},
                "ligand_name": {"existing": "LIG-OLD", "incoming": "LIG-NEW"},
            },
            "existing": [dict(EXISTING)],
        }
    ]


def test_headers_follow_content_fields():
    headers = curation_headers()
    assert headers[0] == "id"
    assert headers[-1] == "action"
    # every content field appears exactly once, in a stable order
    for f in cc.CONTENT_FIELDS:
        assert headers.count(f) == 1
    assert _headers_in(build_curation_xlsx(_conflict_payload())) == headers


def test_untouched_sheet_round_trips_to_defaults():
    data = build_curation_xlsx(_conflict_payload())
    parsed = parse_curation_xlsx(data)

    # incoming row defaults to CREATE, existing duplicate to KEEP
    assert parsed.field_actions[KEY] == "CREATE"
    assert parsed.field_actions[(*KEY, 42)] == "KEEP"
    # merge row is blank until the user picks something
    assert ("MERGE", KEY) not in parsed.field_actions


def test_existing_row_is_keyed_by_incoming_identity_not_its_own():
    """The existing row carries OLD-1 in compound_code; the conflict is exactly
    that it differs from the incoming NEW-1. The decision must still be found
    under the incoming identity, which is what resolve_curation looks up."""
    data = build_curation_xlsx(_conflict_payload())
    parsed = parse_curation_xlsx(data)

    assert (*KEY, 42) in parsed.field_actions
    assert ("AAA", "CCO", "OLD-1", 42) not in parsed.field_actions


def test_delete_decision_reaches_the_resolver():
    payload = _conflict_payload()
    data = _edit(build_curation_xlsx(payload), 42, action="DELETE")
    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert plan.ok
    assert [(a.op, a.existing_id) for a in plan.actions] == [("create", None)]


def test_keep_decision_reaches_the_resolver():
    payload = _conflict_payload()
    plan = resolve_curation(payload, parse_curation_xlsx(build_curation_xlsx(payload)))

    assert plan.ok
    assert [(a.op, a.existing_id) for a in plan.actions] == [("reuse", 42)]


def test_merge_row_is_read_back():
    payload = _conflict_payload()
    data = _edit(
        build_curation_xlsx(payload),
        "Merge",
        action="MERGE",
        compound_code="OLD-1",
        ligand_name="LIG-NEW",
    )
    parsed = parse_curation_xlsx(data)

    merged = parsed.field_actions[("MERGE", KEY)]
    # compound_code is a conflict field like any other and must survive the trip
    assert merged["compound_code"] == "OLD-1"
    assert merged["ligand_name"] == "LIG-NEW"


def test_merge_decision_reaches_the_resolver():
    payload = _conflict_payload()
    data = _edit(
        build_curation_xlsx(payload),
        "Merge",
        action="MERGE",
        compound_code="OLD-1",
        ligand_name="LIG-NEW",
    )
    data = _edit(data, 42, action="MERGE")
    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert plan.ok
    (action,) = plan.actions
    assert action.op == "update"
    assert action.existing_id == 42
    assert action.incoming["compound_code"] == "OLD-1"
    assert action.incoming["ligand_name"] == "LIG-NEW"


def test_blank_merge_cell_fails_rather_than_guessing():
    """A merge-row cell the sheet asked the user to decide must not fall back to
    the incoming value -- that would silently pick a winner on their behalf."""
    payload = _conflict_payload()
    # user chose MERGE but never picked a compound_code from the dropdown
    data = _edit(
        build_curation_xlsx(payload), "Merge", action="MERGE", ligand_name="LIG-OLD"
    )
    data = _edit(data, 42, action="MERGE")
    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert not plan.ok
    assert not plan.actions
    assert "compound_code" in plan.unresolved[0].reason


def test_bad_action_is_rejected():
    data = _edit(build_curation_xlsx(_conflict_payload()), 42, action="DESTROY")
    with pytest.raises(CurationFormatError):
        parse_curation_xlsx(data)


def test_layout_adapts_to_an_added_content_field(monkeypatch):
    """The whole point of the dynamic layout: adding a field to CONTENT_FIELDS
    widens the sheet and round-trips, with no column numbers to update."""
    extra = "assay_comment"
    content = cc.CONTENT_FIELDS + (extra,)
    monkeypatch.setattr(cc, "CONTENT_FIELDS", content)
    monkeypatch.setattr(
        cc,
        "CONFLICT_FIELDS",
        tuple(f for f in content if f not in cc.LOCKED_IDENTITY_FIELDS),
    )

    payload = _conflict_payload()
    payload[0]["incoming"][extra] = "clean"
    payload[0]["existing"][0][extra] = "precipitate"
    payload[0]["conflicts"][extra] = {
        "existing": "precipitate",
        "incoming": "clean",
    }

    data = build_curation_xlsx(payload)
    assert extra in _headers_in(data)

    data = _edit(data, "Merge", action="MERGE", **{extra: "precipitate"})
    parsed = parse_curation_xlsx(data)

    # decisions on the pre-existing columns still land correctly...
    assert parsed.field_actions[KEY] == "CREATE"
    assert parsed.field_actions[(*KEY, 42)] == "KEEP"
    # ...and the new column is read back rather than silently dropped
    assert parsed.field_actions[("MERGE", KEY)][extra] == "precipitate"


def test_group_without_incoming_row_is_rejected():
    data = build_curation_xlsx(_conflict_payload())
    wb = load_workbook(io.BytesIO(data))
    ws = wb[SHEET]
    for row in ws.iter_rows():
        if str(row[0].value).strip().lower() == "incoming":
            ws.delete_rows(row[0].row)
            break
    buf = io.BytesIO()
    wb.save(buf)

    with pytest.raises(CurationFormatError):
        parse_curation_xlsx(buf.getvalue())


def _ambiguous_payload(ids=(39, 40)):
    """Two existing compounds, same structure, different compound_code -- the
    shape of the blocks in a real curation sheet."""
    return [
        {
            "status": "ambiguous",
            "inchi_key": "AAA",
            "incoming": dict(INCOMING),
            "conflicts": {},
            "existing": [
                {
                    "id": ex_id,
                    "inchi_key": "AAA",
                    "smiles": "CCO",
                    "compound_code": f"OLD-{ex_id}",
                    "ligand_name": "LIG-NEW",
                    "modeled_smiles_canon": "CCO",
                }
                for ex_id in ids
            ],
        }
    ]


def test_merge_with_all_actions_deleted_supersedes_the_old_rows():
    """Case 1: DELETE, DELETE, DELETE, MERGE with the merge row filled in."""
    payload = _ambiguous_payload()
    data = build_curation_xlsx(payload)
    for ex_id in (39, 40):
        data = _edit(data, ex_id, action="DELETE")
    data = _edit(data, "Incoming", action="DELETE")
    data = _edit(data, "Merge", action="MERGE", compound_code="CHOSEN")

    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert plan.ok
    (action,) = plan.actions
    assert action.incoming["compound_code"] == "CHOSEN"
    # one compound survives, carrying the merge row's values; the other retires
    assert action.existing_id == 39
    assert action.superseded_ids == [40]


def test_merge_overrides_default_keep_and_create_actions():
    """Case 2: only MERGE selected, the other rows left at their defaults."""
    payload = _ambiguous_payload(ids=(46, 47))
    data = build_curation_xlsx(payload)
    data = _edit(data, "Merge", action="MERGE", compound_code="CHOSEN")

    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert plan.ok
    # exactly one decision -- not a reuse from the defaulted KEEP *and* a create
    (action,) = plan.actions
    assert action.op == "update"
    assert action.incoming["compound_code"] == "CHOSEN"
    assert action.existing_id == 46
    assert action.superseded_ids == [47]


def test_merge_with_unfilled_conflict_fields_is_unresolved():
    """Case 3: MERGE selected but the conflicting cells left blank. This must
    fail the upload rather than quietly picking a winner."""
    payload = _ambiguous_payload(ids=(53, 54))
    # a second conflicting field, so both are left for the user to decide
    payload[0]["existing"][0]["modeled_smiles_soakdb"] = "CCO"
    payload[0]["existing"][1]["modeled_smiles_soakdb"] = "OCC"
    payload[0]["incoming"]["modeled_smiles_soakdb"] = "C-C-O"

    data = build_curation_xlsx(payload)
    data = _edit(data, "Merge", action="MERGE")  # no values chosen

    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert not plan.ok
    assert not plan.actions
    reason = plan.unresolved[0].reason
    assert "compound_code" in reason and "modeled_smiles_soakdb" in reason


def test_merge_with_only_some_fields_filled_is_unresolved():
    payload = _ambiguous_payload(ids=(53, 54))
    payload[0]["existing"][0]["modeled_smiles_soakdb"] = "CCO"
    payload[0]["existing"][1]["modeled_smiles_soakdb"] = "OCC"
    payload[0]["incoming"]["modeled_smiles_soakdb"] = "C-C-O"

    data = build_curation_xlsx(payload)
    data = _edit(data, "Merge", action="MERGE", compound_code="CHOSEN")

    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert not plan.ok
    assert "modeled_smiles_soakdb" in plan.unresolved[0].reason
    assert "compound_code" not in plan.unresolved[0].reason


def test_required_merge_cells_are_marked_and_editable():
    """The cells the user must fill are unlocked and flagged red; the rest are
    pre-filled and locked."""
    payload = _ambiguous_payload()
    ws = load_workbook(io.BytesIO(build_curation_xlsx(payload)))[SHEET]

    cols, merge_row = {}, None
    for row in ws.iter_rows():
        if row[0].value == "id":
            cols = {c.value: c.column for c in row if c.value}
        elif str(row[0].value).strip().lower() == "merge":
            merge_row = row[0].row

    required = ws.cell(merge_row, cols["compound_code"])
    settled = ws.cell(merge_row, cols["ligand_name"])

    assert required.protection.locked is False
    assert required.fill.start_color.rgb == "FFFFC7CE"  # light red
    assert required.value in (None, "")
    assert settled.protection.locked is True
    assert settled.value == "LIG-NEW"


def test_keeping_two_existing_rows_is_unresolved():
    payload = _ambiguous_payload()
    # both default to KEEP, which cannot be satisfied by a single compound
    plan = resolve_curation(payload, parse_curation_xlsx(build_curation_xlsx(payload)))

    assert not plan.ok
    assert "KEEP" in plan.unresolved[0].reason


def test_deleting_every_row_creates_incoming_and_supersedes():
    payload = _ambiguous_payload()
    data = build_curation_xlsx(payload)
    for ex_id in (39, 40):
        data = _edit(data, ex_id, action="DELETE")

    plan = resolve_curation(payload, parse_curation_xlsx(data))

    assert plan.ok
    (action,) = plan.actions
    assert action.op == "create"
    assert action.existing_id is None
    assert sorted(action.superseded_ids) == [39, 40]


def test_multiple_groups_are_kept_separate():
    payload = _conflict_payload()
    payload.append(
        {
            "status": "conflict",
            "inchi_key": "BBB",
            "incoming": {"smiles": "CCC", "compound_code": "NEW-2"},
            "conflicts": {"compound_code": {"existing": "OLD-2", "incoming": "NEW-2"}},
            "existing": [
                {
                    "id": 99,
                    "inchi_key": "BBB",
                    "smiles": "CCC",
                    "compound_code": "OLD-2",
                }
            ],
        }
    )
    data = _edit(build_curation_xlsx(payload), 99, action="DELETE")
    parsed = parse_curation_xlsx(data)

    assert parsed.field_actions[(*KEY, 42)] == "KEEP"
    assert parsed.field_actions[("BBB", "CCC", "NEW-2", 99)] == "DELETE"
