"""Reconcile incoming (upload) compounds against existing ``Compound`` rows.

This is the single source of truth for "does this incoming compound already
exist, and if so is there a conflict". It is deliberately shared by:

* the pre-upload validation endpoint (``UploadExperimentValidateView``) - run
  read-only, *advisory*, so the uploader can warn the user before uploading; and
* the target loader - run for real, *authoritative*, because the database can
  change between validation and upload.

Both callers MUST go through :func:`reconcile_compounds` so the pre-flight check
can never drift from the real load.

Match key: ``(project, inchi_key)`` where ``inchi_key`` is the stereo-preserving
key computed *here* from the incoming ``smiles`` with the same call the loader
uses (``MolFromSmiles`` -> ``MolToInchiKey``, no ``RemoveStereochemistry``).
Callers pass raw compound data (smiles + fields), never a precomputed key - the
key algorithm lives in exactly one place so two codebases can't disagree.
``compound_code`` is deliberately NOT part of the key.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem

from viewer.models import Compound

logger = logging.getLogger(__name__)

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


class MatchStatus(str, Enum):
    CREATE = "create"  # no existing match -> auto-create, no user action
    REUSE = "reuse"  # exactly one match, no field conflict -> auto-reuse
    CONFLICT = "conflict"  # exactly one match, >=1 field differs -> needs user
    AMBIGUOUS = "ambiguous"  # >=2 matches -> needs user


def inchi_key_for_smiles(smiles: str | None) -> str:
    """Stereo-preserving InChIKey, matching the loader exactly.

    Returns "" for missing/unparseable smiles (the loader stores "" too). An
    empty key never matches, so such compounds always classify as CREATE.
    """
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return Chem.inchi.MolToInchiKey(mol)


def _present(value) -> bool:
    """True if the value is real content (not NULL / not blank)."""
    if value is None:
        return False
    if isinstance(value, str) and value.strip() == "":
        return False
    return True


@dataclass
class ExistingCompound:
    id: int
    values: dict  # CONTENT_FIELDS -> stored value


@dataclass
class CompoundMatch:
    inchi_key: str
    incoming: dict  # CONTENT_FIELDS -> supplied value
    status: MatchStatus
    existing: list  # list[ExistingCompound]
    # field -> {"existing": ..., "incoming": ...}, only for a CONFLICT (1 match)
    conflicts: dict = field(default_factory=dict)

    def as_payload(self) -> dict:
        return {
            "inchi_key": self.inchi_key,
            "status": self.status.value,
            "incoming": self.incoming,
            "existing": [{"id": e.id, **e.values} for e in self.existing],
            "conflicts": self.conflicts,
        }


@dataclass
class ReconciliationResult:
    matches: list  # list[CompoundMatch]

    @property
    def needs_curation(self) -> bool:
        return any(
            m.status in (MatchStatus.CONFLICT, MatchStatus.AMBIGUOUS)
            for m in self.matches
        )

    def curation_payload(self) -> list:
        """The subset needing user action, as plain dicts (JSON/spreadsheet)."""
        return [
            m.as_payload()
            for m in self.matches
            if m.status in (MatchStatus.CONFLICT, MatchStatus.AMBIGUOUS)
        ]


def _field_conflicts(existing: Compound, incoming: dict) -> dict:
    """Fields where existing and incoming both have (differing) real content.

    NULL/blank on either side is not a conflict - there is nothing to reconcile,
    the present value simply wins.
    """
    out = {}
    for f in CONTENT_FIELDS:
        ev = getattr(existing, f)
        iv = incoming.get(f)
        if _present(ev) and _present(iv) and ev != iv:
            out[f] = {"existing": ev, "incoming": iv}
    return out


def reconcile_compounds(project, incoming_compounds) -> ReconciliationResult:
    """Classify each incoming compound against existing rows in ``project``.

    ``project``: a ``Project`` instance, or ``None`` when the project does not
    exist yet (first upload) - then everything is CREATE.
    ``incoming_compounds``: iterable of dicts carrying at least ``smiles`` plus
    any of :data:`CONTENT_FIELDS`.
    """
    prepared = [
        (inc, inchi_key_for_smiles(inc.get("smiles"))) for inc in incoming_compounds
    ]

    # One query for every existing compound that could match, keyed by inchi_key.
    keys = {key for _, key in prepared if key}
    existing_by_key: dict = {}
    if project is not None and keys:
        for cmpd in Compound.objects.filter(project=project, inchi_key__in=keys):
            existing_by_key.setdefault(cmpd.inchi_key, []).append(cmpd)

    matches = []
    for inc, key in prepared:
        incoming_vals = {f: inc.get(f) for f in CONTENT_FIELDS}
        existing = existing_by_key.get(key, []) if key else []
        existing_recs = [
            ExistingCompound(id=c.pk, values={f: getattr(c, f) for f in CONTENT_FIELDS})
            for c in existing
        ]

        conflicts: dict = {}
        if not existing:
            status = MatchStatus.CREATE
        elif len(existing) == 1:
            conflicts = _field_conflicts(existing[0], incoming_vals)
            status = MatchStatus.CONFLICT if conflicts else MatchStatus.REUSE
        else:
            status = MatchStatus.AMBIGUOUS

        matches.append(
            CompoundMatch(
                inchi_key=key,
                incoming=incoming_vals,
                status=status,
                existing=existing_recs,
                conflicts=conflicts,
            )
        )

    return ReconciliationResult(matches=matches)
