"""Pure (no-DB) helpers for the compound-curation spreadsheet.

This module deliberately imports nothing from Django/`viewer.models` so it can be
unit-tested in isolation. It owns three things:

* the on-the-wire vocabulary shared by the uploader and the backend
  (:data:`CONTENT_FIELDS`, :func:`curation_headers`, the stable row id);
* generating the curation ``.xlsx`` a user fills in
  (:func:`build_curation_xlsx`) and parsing it back
  (:func:`parse_curation_xlsx`); and
* turning a user's decisions plus a *fresh* reconciliation into a concrete plan
  of per-compound actions (:func:`resolve_curation`).

The workbook holds one visible sheet of per-compound groups - existing
duplicates, the incoming row, then a merge row - plus a hidden metadata sheet.
Its column layout is derived from :data:`CONTENT_FIELDS` at both ends, so adding
a field widens the sheet without any column numbers needing to change.

The spreadsheet is the only artefact that crosses the trust boundary, so parsing
reads *only* the decision cells and the identity columns the sheet locks - and
the loader always re-derives the reconciliation itself before applying a plan.
"""

import hashlib
import io
import json
from dataclasses import dataclass, field

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Protection, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

# Compound fields populated from an upload (meta_aligner.yaml), excluding the
# identity key (inchi_key) and the project FK. A difference in any of these
# between an incoming compound and its single existing match is a conflict.
CONTENT_FIELDS = (
    "smiles",
    "compound_code",
    "ligand_name",
    "modeled_smiles_soakdb",
    "modeled_smiles_canon",
    "soaked_smiles_soakdb",
    "soaked_smiles_canon",
)

# Fields that identify a compound (locked in identity columns, can't be changed)
# Note: inchi_key and smiles are pure identity; compound_code can conflict
LOCKED_IDENTITY_FIELDS = ("inchi_key", "smiles")

# Fields that can have conflicts (all CONTENT_FIELDS that aren't pure identity)
CONFLICT_FIELDS = tuple(f for f in CONTENT_FIELDS if f not in LOCKED_IDENTITY_FIELDS)

# Column headings that are not compound fields.
ROW_ID_HEADER = "id"
ACTION_HEADER = "action"

# Written into the id column to mark the two synthetic rows of a group. Matched
# case-insensitively on the way back in.
ROW_MARKER_INCOMING = "Incoming"
ROW_MARKER_MERGE = "Merge"


def curation_headers() -> list[str]:
    """Canonical column order for one compound-group subtable.

    Derived from CONTENT_FIELDS, so adding a field there widens the sheet with
    no other change: the writer lays cells out by looking a name up in this
    list, and the parser rebuilds the same map from the header row it reads
    back. Neither side hardcodes a column number.
    """
    return [ROW_ID_HEADER, *LOCKED_IDENTITY_FIELDS, *CONFLICT_FIELDS, ACTION_HEADER]


def _column_map(headers) -> dict:
    """Map column name -> 0-based index, for indexing a ``values_only`` row."""
    return {name: idx for idx, name in enumerate(headers) if name}


def merge_candidate_values(name, existing, incoming) -> list:
    """Distinct non-empty values on offer for ``name`` in a group's merge row.

    Order is existing-rows-first, then incoming, so the dropdown reads in the
    same order as the rows above it.
    """
    values: list = []
    for source in (*(existing or []), incoming or {}):
        value = (source or {}).get(name)
        if value not in (None, "") and value not in values:
            values.append(value)
    return values


def fields_needing_merge_choice(existing, incoming) -> list:
    """Conflict fields the merge row leaves for the user to decide.

    A field only needs a decision when two or more distinct non-empty values
    are on offer; with one (or none) the writer pre-fills and locks the cell.
    The writer uses this to decide which cells get a dropdown, and
    :func:`resolve_curation` uses it to check the user actually filled them in -
    so the sheet and the validation cannot disagree about what "complete" means.
    """
    return [
        name
        for name in CONFLICT_FIELDS
        if len(merge_candidate_values(name, existing, incoming)) > 1
    ]


def _cell(row, cols, name):
    """Value of a named column in a ``values_only`` row tuple, or None.

    Tolerates both a short row (openpyxl trims trailing empties) and a column
    the sheet simply does not carry.
    """
    idx = cols.get(name)
    if idx is None or idx >= len(row):
        return None
    return row[idx]


STATUS_CONFLICT = "conflict"
STATUS_AMBIGUOUS = "ambiguous"
# A conflict auto_merge_non_conflicting_compounds settled on its own: the user
# is never asked about it, but it still carries a decision the loader must apply
# (reuse the matched row and write the merged content onto it).
STATUS_AUTO_MERGED = "auto_merged"

# Statuses that still need a human decision, and so appear in the spreadsheet.
CURATION_STATUSES = (STATUS_CONFLICT, STATUS_AMBIGUOUS)

META_SHEET = "_curation"

_LOCKED_FILL = PatternFill(
    start_color="FFEFEFEF", end_color="FFEFEFEF", fill_type="solid"
)
_EDIT_FILL = PatternFill(
    start_color="FFFFF3C4", end_color="FFFFF3C4", fill_type="solid"
)
# A merge-row cell the user MUST fill in: light red, so an unresolved conflict
# is visible at a glance rather than only surfacing as an upload failure.
_REQUIRED_FILL = PatternFill(
    start_color="FFFFC7CE", end_color="FFFFC7CE", fill_type="solid"
)
_HEADER_FONT = Font(bold=True)
_GREY_FONT = Font(name="Calibri", size=11, color="FF999999", bold=False, italic=False)
_WRAP = Alignment(wrap_text=True, vertical="top")


class CurationFormatError(Exception):
    """The uploaded curation spreadsheet is malformed or self-contradictory."""


def compound_row_id(inchi_key, incoming) -> str:
    """Stable opaque id for one incoming compound.

    Keyed by (inchi_key, incoming smiles, incoming compound_code) so a re-run
    against the same upload reproduces the same ids, letting the loader match a
    user's decision to a freshly-derived conflict without trusting any editable
    cell.
    """
    basis = "\x00".join(
        [
            inchi_key or "",
            (incoming or {}).get("smiles") or "",
            (incoming or {}).get("compound_code") or "",
        ]
    )
    return hashlib.sha1(basis.encode("utf-8")).hexdigest()[:16]


def curation_payload_hash(payload) -> str:
    """Content hash of a curation payload, for provenance / sanity checks."""
    blob = json.dumps(payload, sort_keys=True, default=str, ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def _lock(cell):
    cell.protection = Protection(locked=True)
    cell.fill = _LOCKED_FILL
    cell.alignment = _WRAP


def _editable(cell):
    cell.protection = Protection(locked=False)
    cell.fill = _EDIT_FILL


def auto_merge_non_conflicting_compounds(payload) -> tuple[list, dict]:
    """Settle the conflicts that need no human input, in place.

    Criteria:
    - NULL/empty values don't create conflicts (existing value wins)
    - For fields like description/comments, concatenate instead of flagging conflict
    - Only genuine conflicts (both sides have different non-empty values) are kept

    Returns (payload, stats). The payload has the SAME length and order as the
    input: a conflict this resolves is rewritten to STATUS_AUTO_MERGED carrying
    the merged values, not removed. It has to stay, because it still describes a
    decision - reuse the matched row - that the loader must apply; dropping it
    left the loader with no action and it created a duplicate compound instead.
    Use :func:`needs_curation` to get just the entries a user must act on.
    """
    CONCATENABLE_FIELDS = {"description", "comments"}

    curation_payload = []
    stats: dict = {"auto_merged": 0, "genuine_conflicts": 0, "fields_consolidated": {}}

    for item in payload:
        if item.get("status") != STATUS_CONFLICT:
            curation_payload.append(item)
            continue

        incoming = item.get("incoming", {})
        existing_list = item.get("existing", [])
        conflicts = item.get("conflicts", {})

        # Filter conflicts and build merged values
        genuine_conflicts = {}
        merged_values = {}

        for field_name, conflict_info in conflicts.items():
            existing_val = conflict_info.get("existing", "")
            incoming_val = conflict_info.get("incoming", "")

            # Normalize to string and strip for comparison
            existing_str = str(existing_val).strip() if existing_val else ""
            incoming_str = str(incoming_val).strip() if incoming_val else ""

            # For concatenable fields, combine values
            if field_name in CONCATENABLE_FIELDS:
                parts = []
                if existing_str:
                    parts.append(existing_str)
                if incoming_str:
                    parts.append(incoming_str)
                if parts:
                    merged_values[field_name] = " | ".join(parts)
                    stats["fields_consolidated"][field_name] = (
                        stats["fields_consolidated"].get(field_name, 0) + 1
                    )
                continue

            # For other fields: NULL/empty doesn't create conflict
            if not existing_str and incoming_str:
                # Incoming has value, existing is empty -> use incoming
                merged_values[field_name] = incoming_val
                continue
            if existing_str and not incoming_str:
                # Existing has value, incoming is empty -> use existing
                merged_values[field_name] = existing_val
                continue
            if not existing_str and not incoming_str:
                # Both empty -> no conflict
                continue

            # Both sides have values and they differ -> genuine conflict
            if existing_str != incoming_str:
                genuine_conflicts[field_name] = conflict_info

        # If no genuine conflicts remain, auto-merge
        if not genuine_conflicts:
            stats["auto_merged"] += 1
            # Update existing compound with merged values
            if existing_list:
                existing_list[0].update(merged_values)
            # Update incoming with merged values for consistency
            merged_incoming = dict(incoming)
            merged_incoming.update(merged_values)
            item_copy = dict(item)
            item_copy["incoming"] = merged_incoming
            item_copy["status"] = STATUS_AUTO_MERGED
            # Kept, not dropped: resolve_curation turns this into the "reuse the
            # matched row" action. build_curation_xlsx ignores the status, so it
            # still never reaches the spreadsheet.
            curation_payload.append(item_copy)
        else:
            # Keep for manual curation
            stats["genuine_conflicts"] += 1
            item_copy = dict(item)
            item_copy["conflicts"] = genuine_conflicts
            curation_payload.append(item_copy)

    return curation_payload, stats


def needs_curation(payload) -> list:
    """The payload entries that still require a human decision."""
    return [m for m in payload if m.get("status") in CURATION_STATUSES]


def build_curation_xlsx(curation_payload, *, target_name=None) -> bytes:
    """Render conflicts/duplicates into a hierarchical compound-group workbook.

    Each incoming compound that has conflicts or duplicates gets its own subtable:
    - Existing duplicates (with conflicting fields)
    - Incoming compound row
    - Merge row (where user specifies which values win, only for conflicting fields)
    """
    conflicts = [m for m in curation_payload if m.get("status") == STATUS_CONFLICT]
    duplicates = [m for m in curation_payload if m.get("status") == STATUS_AMBIGUOUS]

    wb = Workbook()

    # Single sheet: Incoming compounds with conflicts and duplicates in compound groups
    _build_incoming_molecules_sheet(wb.active, conflicts, duplicates)
    wb.active.title = "Incoming compounds"

    # Hidden metadata sheet
    _build_meta_sheet(wb.create_sheet(META_SHEET), curation_payload, target_name)

    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def _protect(ws):
    ws.protection.sheet = True
    ws.protection.selectLockedCells = (
        True  # Allow selecting locked cells with validation dropdowns
    )
    ws.protection.enable()


def _apply_group_border(ws, start_row, end_row, inner_border, outer_border):
    """Apply borders to group rows for visual separation."""
    for row_num in range(start_row, end_row + 1):
        for col in range(1, ws.max_column + 1):
            cell = ws.cell(row_num, col)

            top = outer_border if row_num == start_row else inner_border
            bottom = outer_border if row_num == end_row else inner_border
            left = outer_border if col == 1 else inner_border
            right = outer_border if col == ws.max_column else inner_border

            cell.border = Border(top=top, bottom=bottom, left=left, right=right)


def _build_incoming_molecules_sheet(ws, conflicts, duplicates):
    """Build sheet with compound groups: existing dups, incoming, merge row per group.

    Merge row is populated only with NON-conflicting fields (fields that are identical
    across existing and incoming compounds). Conflicting fields are left empty for users
    to decide which value to keep when selecting the MERGE action.
    """
    ws.title = "incoming molecules"

    # Merge conflicts and duplicates into compound groups by incoming compound
    groups = (
        {}
    )  # (inchi_key, incoming_smiles, incoming_code) -> {existing: [...], conflicts: {...}}
    for m in conflicts:
        incoming = m.get("incoming", {})
        key = (
            m.get("inchi_key", ""),
            incoming.get("smiles", ""),
            incoming.get("compound_code", ""),
        )
        if key not in groups:
            groups[key] = {
                "existing": [],
                "conflicts": m.get("conflicts", {}),
                "incoming": incoming,
            }
        groups[key]["existing"].append(m.get("existing", [{}])[0])

    for m in duplicates:
        incoming = m.get("incoming", {})
        key = (
            m.get("inchi_key", ""),
            incoming.get("smiles", ""),
            incoming.get("compound_code", ""),
        )
        if key not in groups:
            groups[key] = {"existing": [], "conflicts": {}, "incoming": incoming}
        groups[key]["existing"].extend(m.get("existing", []))

    # Every column number below is looked up by name in this one list, so the
    # layout stays consistent however CONTENT_FIELDS grows or is reordered.
    headers = curation_headers()
    field_to_col = {f: headers.index(f) + 1 for f in CONFLICT_FIELDS}
    action_col = headers.index(ACTION_HEADER) + 1

    # Column widths, keyed the same way. SMILES-ish columns get the wide
    # treatment; the two identity columns have hand-tuned widths.
    fixed_widths = {"inchi_key": 28, "smiles": 50}
    col_widths = {1: 15}  # id
    for name in (*LOCKED_IDENTITY_FIELDS, *CONFLICT_FIELDS):
        col_widths[headers.index(name) + 1] = fixed_widths.get(
            name, 40 if "smiles" in name else 18
        )

    row_num = 1
    action_dv = DataValidation(
        type="list",
        formula1='"DELETE,KEEP,CREATE,MERGE"',
        allow_blank=False,
        showDropDown=False,
    )
    ws.add_data_validation(action_dv)

    for key, group_data in groups.items():
        inchi_key, incoming_smiles, _ = key
        existing = group_data["existing"]
        incoming = group_data["incoming"]
        conflicts = group_data["conflicts"]

        # Write header row for this group
        for col, header in enumerate(headers, 1):
            cell = ws.cell(row_num, col, header)
            cell.font = _HEADER_FONT
            cell.fill = _LOCKED_FILL

        row_num += 1

        # Write existing duplicates
        for existing_compound in existing:
            ws.cell(row_num, 1, existing_compound.get("id", ""))
            # Write locked identity fields
            for idx, identity_field in enumerate(LOCKED_IDENTITY_FIELDS, start=2):
                ws.cell(row_num, idx, existing_compound.get(identity_field, ""))
            # Write conflict fields (includes compound_code)
            for name, col_num in field_to_col.items():
                ws.cell(row_num, col_num, existing_compound.get(name, ""))

            action_cell = ws.cell(row_num, action_col, "KEEP" if existing else "DELETE")
            _editable(action_cell)
            action_dv.add(action_cell)

            for col in range(1, action_col):
                _lock(ws.cell(row_num, col))

            row_num += 1

        # Write incoming row (visually distinct with light background)
        incoming_bg = PatternFill(
            start_color="FFF0F0F0", end_color="FFF0F0F0", fill_type="solid"
        )
        id_cell = ws.cell(row_num, 1, ROW_MARKER_INCOMING)
        id_cell.font = Font(bold=True)
        id_cell.fill = incoming_bg

        # Write locked identity fields
        locked_identity_values = [inchi_key, incoming_smiles]
        for idx, value in enumerate(locked_identity_values, start=2):
            ws.cell(row_num, idx, value).fill = incoming_bg
        # Write conflict fields (includes compound_code)
        for name, col_num in field_to_col.items():
            ws.cell(row_num, col_num, incoming.get(name, "")).fill = incoming_bg

        action_cell = ws.cell(row_num, action_col, "CREATE")
        _editable(action_cell)
        action_dv.add(action_cell)

        for col in range(1, action_col):
            _lock(ws.cell(row_num, col))

        row_num += 1

        # Write merge row - only populate NON-conflicting fields
        merge_id = ws.cell(row_num, 1, ROW_MARKER_MERGE)
        merge_id.font = Font(bold=True)  # Keep "Merge" black and bold
        merge_id.protection = Protection(locked=True)
        merge_id.fill = _LOCKED_FILL
        merge_id.alignment = Alignment(vertical="top", wrap_text=True)

        # Write locked identity fields (inchi_key, smiles - never conflict)
        locked_identity_values = [inchi_key, incoming_smiles]
        for idx, value in enumerate(locked_identity_values, start=2):
            cell = ws.cell(row_num, idx, value)
            cell.protection = Protection(locked=True)
            cell.fill = _LOCKED_FILL
            cell.font = _GREY_FONT
            cell.alignment = Alignment(vertical="top", wrap_text=True)

        # Process conflict fields dynamically. Which cells need a decision is
        # decided by the same helper resolve_curation validates against, so the
        # sheet cannot offer a choice the loader then fails to require.
        needs_choice = set(fields_needing_merge_choice(existing, incoming))
        for name in CONFLICT_FIELDS:
            col_num = field_to_col[name]
            candidates = merge_candidate_values(name, existing, incoming)

            if name in needs_choice:
                # Two or more values on offer: the user MUST pick one. Red so it
                # reads as "action required" - leaving it blank fails the upload.
                cell = ws.cell(row_num, col_num, "")
                cell.protection = Protection(locked=False)
                cell.fill = _REQUIRED_FILL
                cell.alignment = Alignment(vertical="top", wrap_text=True)
                options = ",".join(str(v) for v in candidates)
                dv = DataValidation(
                    type="list",
                    formula1=f'"{options}"',
                    allow_blank=True,
                    showDropDown=False,
                )
                ws.add_data_validation(dv)
                dv.add(cell)
            else:
                # Nothing to decide: pre-fill with the single value on offer (or
                # blank if neither side has one) and lock it.
                cell = ws.cell(row_num, col_num, candidates[0] if candidates else "")
                cell.protection = Protection(locked=True)
                cell.fill = _LOCKED_FILL
                cell.font = _GREY_FONT
                cell.alignment = Alignment(vertical="top", wrap_text=True)

        merge_action = ws.cell(row_num, action_col, "")
        _editable(merge_action)
        action_dv.add(merge_action)

        row_num += 1

        # Blank separator rows between groups
        row_num += 2

    # Apply borders around groups
    thin_border = Side(style="thin", color="D3D3D3")
    group_border = Side(style="medium", color="000000")

    current_group_start = None
    for row_num in range(1, ws.max_row + 1):
        row_vals = [ws.cell(row_num, col).value for col in range(1, ws.max_column + 1)]

        # Header row marks start of group
        if row_vals[0] == ROW_ID_HEADER:
            current_group_start = row_num + 1
        # Blank row marks end of group
        elif not row_vals[0] or all(v is None or v == "" for v in row_vals):
            if current_group_start is not None:
                _apply_group_border(
                    ws, current_group_start, row_num - 1, thin_border, group_border
                )
                current_group_start = None

    # Apply borders to last group if any
    if current_group_start is not None:
        _apply_group_border(
            ws, current_group_start, ws.max_row, thin_border, group_border
        )

    # Apply dynamic column widths
    for col_num, width in col_widths.items():
        ws.column_dimensions[get_column_letter(col_num)].width = width
    ws.column_dimensions[get_column_letter(action_col)].width = 12

    _protect(ws)


def _build_meta_sheet(ws, curation_payload, target_name):
    ws.append(["key", "value"])
    ws.append(["payload_hash", curation_payload_hash(curation_payload)])
    ws.append(["target_name", target_name or ""])
    ws.append(
        [
            "do_not_edit",
            "This sheet lets the server identify the conflict set. Do not edit.",
        ]
    )
    ws.sheet_state = "hidden"
    _protect(ws)


@dataclass
class ParsedCuration:
    payload_hash: str
    # keyed as _parse_incoming_molecules_sheet documents: the incoming identity
    # triple, that triple plus an existing row id, or ("MERGE", triple).
    field_actions: dict = field(default_factory=dict)


def parse_curation_xlsx(data: bytes, match_payloads=None) -> ParsedCuration:
    """Read the new hierarchical curation workbook into decisions.

    Parses compound groups from "Incoming compounds" sheet:
    - Existing duplicates with DELETE/KEEP actions
    - Incoming compound row with CREATE action
    - Merge row where user specifies merged values

    Returns field_actions keyed by compound identity (inchi_key, smiles, code).
    """
    wb = load_workbook(io.BytesIO(data), data_only=True)

    payload_hash = ""
    if META_SHEET in wb.sheetnames:
        for row in wb[META_SHEET].iter_rows(min_row=2, values_only=True):
            if row and row[0] == "payload_hash":
                payload_hash = row[1] or ""

    # Parse the hierarchical "Incoming compounds" sheet
    field_actions = {}
    if "Incoming compounds" in wb.sheetnames:
        field_actions = _parse_incoming_molecules_sheet(
            wb["Incoming compounds"], match_payloads
        )

    return ParsedCuration(payload_hash=payload_hash, field_actions=field_actions)


VALID_ROW_ACTIONS = ("DELETE", "KEEP", "CREATE", "MERGE")


def _flush_curation_group(rows, cols, actions) -> None:
    """Turn one compound group's rows into entries in ``actions``.

    ``rows`` is [(sheet_row_number, values_only tuple)] for a single group, in
    sheet order; ``cols`` maps column name -> index for that group.
    """
    if not rows:
        return

    # The whole group is keyed by the *incoming* compound's identity. The
    # existing rows carry their own values in those columns - differing from
    # the incoming ones is precisely what makes them a conflict - so they
    # cannot be used to key the group, and the merge row's compound_code is a
    # user-editable choice rather than an identity.
    key = None
    for _, row in rows:
        if str(row[0]).strip().lower() == ROW_MARKER_INCOMING.lower():
            key = (
                _cell(row, cols, "inchi_key"),
                _cell(row, cols, "smiles"),
                _cell(row, cols, "compound_code"),
            )
            break

    if key is None:
        raise CurationFormatError(
            f"Row {rows[0][0]}: compound group has no '{ROW_MARKER_INCOMING}' row; "
            "the curation file appears to be corrupt."
        )

    for row_num, row in rows:
        marker = str(row[0]).strip().lower()
        action = (_cell(row, cols, ACTION_HEADER) or "").strip()

        if marker == ROW_MARKER_INCOMING.lower():
            if action:
                actions[key] = action

        elif marker == ROW_MARKER_MERGE.lower():
            if action.upper() == "MERGE":
                # Read back every conflict field the sheet actually carries, so
                # a decision on a newly-added CONTENT_FIELD is not silently
                # dropped. Blank cells (a dropdown the user never picked) are
                # omitted so the incoming value stands rather than being
                # clobbered with None.
                actions[("MERGE", key)] = {
                    f: _cell(row, cols, f)
                    for f in CONFLICT_FIELDS
                    if f in cols and _cell(row, cols, f) not in (None, "")
                }

        else:
            # An existing duplicate row; its id column holds the database id.
            if not action:
                continue
            if action not in VALID_ROW_ACTIONS:
                raise CurationFormatError(
                    f"Row {row_num}: '{action}' is not a valid action; "
                    f"expected {', '.join(VALID_ROW_ACTIONS)}"
                )
            actions[(*key, row[0])] = action


def _parse_incoming_molecules_sheet(  # pylint: disable=unused-argument
    ws, match_payloads=None
) -> dict:
    """Parse compound groups from the "Incoming compounds" sheet.

    Column positions come from each group's own header row rather than being
    assumed, so a sheet built from a wider or reordered CONTENT_FIELDS reads
    back correctly without changing anything here.

    Returns an actions dict keyed by the group's incoming identity triple
    (inchi_key, smiles, compound_code): the bare triple for the incoming row,
    the triple plus the existing row's database id for each existing duplicate,
    and ("MERGE", triple) for the merge row's chosen values.

    ``match_payloads`` is accepted for symmetry with the other sheet parsers but
    is not consulted: every identity column here is locked, and the loader
    re-derives the reconciliation from the database before applying anything.
    """
    actions: dict = {}
    cols = _column_map(curation_headers())
    group: list = []

    for row_num, row in enumerate(ws.iter_rows(min_row=1, values_only=True), start=1):
        if not row or row[0] in (None, ""):
            continue

        if row[0] == ROW_ID_HEADER:
            # A header row both closes the previous group and describes the
            # columns of the next one.
            _flush_curation_group(group, cols, actions)
            cols = _column_map(row)
            group = []
            continue

        group.append((row_num, row))

    _flush_curation_group(group, cols, actions)

    return actions


@dataclass
class ResolvedAction:
    index: int  # position in the incoming-compound list
    op: str  # "create" | "reuse" | "update"
    existing_id: int | None
    incoming: dict
    # Existing compounds this decision retires: their references are relinked
    # to the surviving compound and the rows themselves deleted.
    superseded_ids: list = field(default_factory=list)


@dataclass
class Unresolved:
    index: int
    inchi_key: str
    reason: str


@dataclass
class ResolutionPlan:
    actions: list  # list[ResolvedAction]
    unresolved: list  # list[Unresolved]

    @property
    def ok(self) -> bool:
        return not self.unresolved


def _resolve_group(i, m, parsed, actions, unresolved) -> None:
    """Turn one curated compound group into at most one action.

    A group yields exactly one decision, because it describes a single incoming
    compound. MERGE on the merge row wins over whatever the per-row actions say:
    it is the only decision that can express "fold all of these into one", so
    the KEEP/DELETE cells left at their defaults must not contradict it.
    """
    inchi_key = m.get("inchi_key", "")
    incoming = m.get("incoming", {})
    existing = m.get("existing", []) or []
    key = (inchi_key, incoming.get("smiles", ""), incoming.get("compound_code", ""))

    existing_ids = [e["id"] for e in existing if e.get("id") is not None]

    merge_values = parsed.field_actions.get(("MERGE", key))
    if merge_values is not None:
        # The user chose MERGE. Every field the sheet asked them to decide must
        # actually carry a value - a blank one means the conflict is unresolved,
        # and guessing would silently pick a winner on their behalf.
        missing = [
            name
            for name in fields_needing_merge_choice(existing, incoming)
            if merge_values.get(name) in (None, "")
        ]
        if missing:
            unresolved.append(
                Unresolved(
                    i,
                    inchi_key,
                    "MERGE was selected but no value was chosen for: "
                    + ", ".join(missing),
                )
            )
            return

        merged_incoming = dict(incoming)
        merged_incoming.update(merge_values)
        # The lowest existing id survives and takes the merged content; the rest
        # are relinked onto it and deleted.
        survivor = min(existing_ids) if existing_ids else None
        superseded = [e for e in existing_ids if e != survivor]
        actions.append(
            ResolvedAction(i, "update", survivor, merged_incoming, superseded)
        )
        return

    # No MERGE: each existing row's own KEEP/DELETE decides its fate.
    per_row = {
        ex_id: parsed.field_actions.get((*key, ex_id), "KEEP") for ex_id in existing_ids
    }
    kept = [ex_id for ex_id, act in per_row.items() if act == "KEEP"]
    deleted = [ex_id for ex_id, act in per_row.items() if act == "DELETE"]

    if len(kept) > 1:
        unresolved.append(
            Unresolved(
                i,
                inchi_key,
                f"{len(kept)} existing compounds are marked KEEP; mark all but one "
                "DELETE, or use MERGE to fold them into a single compound",
            )
        )
        return

    if kept:
        # Reuse the kept row; anything explicitly deleted folds into it.
        actions.append(ResolvedAction(i, "reuse", kept[0], incoming, deleted))
        return

    if deleted:
        # Everything retired: the incoming compound is created fresh and the
        # retired rows are relinked onto it.
        actions.append(ResolvedAction(i, "create", None, incoming, deleted))
        return

    unresolved.append(Unresolved(i, inchi_key, "compound has no curation decision"))


def resolve_curation(match_payloads, parsed: ParsedCuration | None) -> ResolutionPlan:
    """Combine a fresh reconciliation with the user's curation decisions.

    Compounds that reconciled cleanly (create/reuse) need no decision. For a
    conflict or an ambiguous multi-match the user's sheet decides, per group:

    - MERGE on the merge row: fold the group into one compound carrying the
      merge row's values. Overrides the per-row actions. Every field the sheet
      offered a choice for must be filled in, or the group is left unresolved.
    - KEEP on exactly one existing row: reuse it; any DELETE rows fold into it.
    - DELETE on all existing rows: create the incoming compound fresh and fold
      the deleted rows into it.

    Anything else is reported through :attr:`ResolutionPlan.unresolved`, which
    the loader turns into an error that rolls the whole upload back.
    """
    actions: list = []
    unresolved: list = []

    for i, m in enumerate(match_payloads):
        status = m.get("status")
        incoming = m.get("incoming", {})
        existing = m.get("existing", [])

        if status == "create":
            actions.append(ResolvedAction(i, "create", None, incoming))

        elif status == "reuse":
            actions.append(ResolvedAction(i, "reuse", existing[0]["id"], incoming))

        elif status == STATUS_AUTO_MERGED:
            # Settled without asking the user: reuse the matched row and write
            # the merged content onto it. `incoming` is already the merged dict.
            existing_id = existing[0]["id"] if existing else None
            actions.append(
                ResolvedAction(
                    i,
                    "update" if existing_id is not None else "create",
                    existing_id,
                    incoming,
                )
            )

        elif status in CURATION_STATUSES:
            if not parsed:
                unresolved.append(
                    Unresolved(
                        i,
                        m.get("inchi_key", ""),
                        f"{status} compound has no curation decision",
                    )
                )
                continue
            _resolve_group(i, m, parsed, actions, unresolved)

        else:
            unresolved.append(
                Unresolved(
                    i,
                    m.get("inchi_key", ""),
                    f"unknown reconciliation status {status!r}",
                )
            )

    return ResolutionPlan(actions=actions, unresolved=unresolved)
