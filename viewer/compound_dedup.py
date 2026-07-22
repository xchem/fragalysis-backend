"""Investigate / collapse duplicate Compound rows.

    ANALYSIS TOOL - NOT PART OF THE APPLICATION.
    DO NOT WRITE ANY CODE THAT REQUIRES THIS MODULE.
    This module and its ``dedup_compounds`` management command exist only to
    investigate duplicate compounds as part of 2300. Actual
    deduplication was decided to require human curation, so nothing here is
    wired into the app. It is safe to delete this file and
    ``viewer/management/commands/dedup_compounds.py`` at any time.

Now that ``Compound.project`` is a scalar FK, duplicates are the rows sharing a
``(project, inchi_key, smiles, compound_code)`` - i.e. same project and an exact
match on structure key, SMILES string and code. Each such group should collapse
to exactly one canonical row (the lowest-pk member), with *every* object that
points at a duplicate - through any FK or M2M - repointed to the survivor.

The grouping key is intentionally conservative: including the raw ``smiles`` and
``compound_code`` means rows only merge when all three agree, so nothing is
over-merged while the exact dedup rules are still being decided.

DRY RUN: :func:`collapse_duplicate_compounds` does NOT modify the database. It
returns a :class:`DedupReport` describing what would happen (how many rows would
be removed, and how many referring objects each relation would move). The actual
repoint/delete implementation lives in :func:`_repoint` with its mutating calls
commented out, ready to enable.

The report also flags rows that share the full ``(project, inchi_key, smiles,
compound_code)`` key but disagree on some *other* content field (ignoring
NULL/blank values) - a heads-up that rows about to collapse together carry
divergent data worth inspecting before committing to any collapse rule.

Referrers are discovered from ``Compound._meta`` at runtime, so nothing is missed
- on the current schema that is SiteObservation, Pose, Result, ActivityPoint,
scoring.CmpdChoice, hypothesis.Vector, the ExperimentCompound / CompoundIdentifier
join rows, the DesignSet.compounds M2M, and the inspirations M2M.
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field

from viewer.models import Compound

# This module inspects the model graph via Django's ``_meta`` API to discover
# every relation pointing at Compound - that's the intended way to do it.
# pylint: disable=protected-access

logger = logging.getLogger(__name__)


@dataclass
class DedupReport:
    total_compounds: int = 0
    groups: int = 0
    duplicate_groups: int = 0
    would_remove: int = 0
    # relation label -> number of referring objects that would move to the keeper
    fk_referrers: dict = field(default_factory=dict)
    m2m_referrers: dict = field(default_factory=dict)
    # forward m2m on Compound whose links would be folded into the keeper
    would_fold: dict = field(default_factory=dict)
    # sample of duplicate groups: (project_id, inchi_key, smiles, code, count)
    sample_groups: list = field(default_factory=list)
    # rows sharing (project, inchi_key, smiles, compound_code) but disagreeing
    # on some *other* content field. field name -> number of conflict groups.
    field_conflict_groups: int = 0
    field_conflict_field_counts: dict = field(default_factory=dict)
    # detail: (project_id, inchi_key, smiles, code, {field: [distinct values]})
    field_conflicts: list = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            "",
            "=== Compound dedup (DRY RUN - nothing written) ===",
            f"  compounds:                                  {self.total_compounds}",
            f"  (project, inchi_key, smiles, code) groups:  {self.groups}",
            f"  duplicate groups:                           {self.duplicate_groups}"
            f"  -> would remove {self.would_remove} rows",
            "  referrers that would repoint to the keeper (FK):",
        ]
        lines += [f"      {k}: {v}" for k, v in sorted(self.fk_referrers.items())]
        lines.append("  referrers via M2M:")
        lines += [f"      {k}: {v}" for k, v in sorted(self.m2m_referrers.items())]
        lines.append("  keeper absorbs the duplicates' own M2M links:")
        lines += [f"      {k}: {v}" for k, v in sorted(self.would_fold.items())]
        lines.append(
            "  same (project, inchi_key, smiles, code) but conflicting content:"
            f"  {self.field_conflict_groups} groups"
        )
        lines += [
            f"      {k}: {v}"
            for k, v in sorted(self.field_conflict_field_counts.items())
        ]
        return "\n".join(lines)


def _through_compound_fk_attname(through) -> str:
    """The attname of the FK on ``through`` that points at Compound."""
    for f in through._meta.get_fields():
        if getattr(f, "related_model", None) is Compound and not f.many_to_many:
            return f.attname
    raise LookupError(f"no Compound FK on {through._meta.label}")


def _content_field_names() -> list:
    """Concrete, non-relational Compound fields excluding the grouping keys.

    These are the 'other content' fields checked for conflicts within a
    (project, inchi_key, smiles, compound_code) group - inchi, descriptions,
    ligand_name, the modeled/soaked smiles, etc. The grouping keys, relations
    (project, current_identifier, inspirations) and the pk are excluded.
    """
    grouping = {"inchi_key", "smiles", "compound_code"}
    return [
        f.attname
        for f in Compound._meta.concrete_fields
        if not f.primary_key and not f.is_relation and f.name not in grouping
    ]


def _is_present(value) -> bool:
    """True if the value is real content (not NULL / not blank)."""
    if value is None:
        return False
    if isinstance(value, str) and value.strip() == "":
        return False
    return True


def _collect_field_conflicts(report: "DedupReport") -> None:
    """Find rows that share (project, inchi_key, smiles, compound_code) but
    disagree on any other content field (ignoring NULL/blank values)."""
    attnames = _content_field_names()
    groups: dict = defaultdict(list)
    for row in (
        Compound.objects.exclude(inchi_key="")
        .exclude(inchi_key__isnull=True)
        .values("project_id", "inchi_key", "smiles", "compound_code", *attnames)
        .iterator()
    ):
        key = (
            row["project_id"],
            row["inchi_key"],
            row["smiles"],
            row["compound_code"],
        )
        groups[key].append(row)

    tally: dict = defaultdict(int)
    for key, members in groups.items():
        if len(members) < 2:
            continue
        conflicts: dict = {}
        for attn in attnames:
            present = {m[attn] for m in members if _is_present(m[attn])}
            if len(present) >= 2:
                conflicts[attn] = sorted(present, key=str)
        if conflicts:
            report.field_conflict_groups += 1
            for attn in conflicts:
                tally[attn] += 1
            report.field_conflicts.append((*key, conflicts))
    report.field_conflict_field_counts = dict(tally)


def collapse_duplicate_compounds() -> DedupReport:
    report = DedupReport()
    report.total_compounds = Compound.objects.count()

    groups: dict = defaultdict(list)
    for cmpd in (
        Compound.objects.exclude(inchi_key="")
        .exclude(inchi_key__isnull=True)
        .only("pk", "project_id", "inchi_key", "smiles", "compound_code")
    ):
        key = (cmpd.project_id, cmpd.inchi_key, cmpd.smiles, cmpd.compound_code)
        groups[key].append(cmpd)
    report.groups = len(groups)

    dup_pks: list = []
    for key, members in groups.items():
        if len(members) < 2:
            continue
        report.duplicate_groups += 1
        project_id, inchi_key, smiles, code = key
        report.sample_groups.append((project_id, inchi_key, smiles, code, len(members)))

        # every member is an exact (project, inchi_key, smiles, code) match, so
        # the whole group collapses to the lowest-pk survivor.
        keeper = min(members, key=lambda c: c.pk)
        for dup in members:
            if dup.pk == keeper.pk:
                continue
            dup_pks.append(dup.pk)
            # --- COLLAPSE (disabled - nothing is written yet) ---
            # _repoint(dup, keeper)
            # dup.delete()
    report.would_remove = len(dup_pks)

    # Tally every referrer of the duplicates, read-only, so the full blast radius
    # is visible before we commit to a collapse strategy.
    for rel in Compound._meta.related_objects:
        label = f"{rel.related_model._meta.label}.{rel.field.name}"
        if rel.many_to_many:
            through = rel.through
            attn = _through_compound_fk_attname(through)
            report.m2m_referrers[label] = through.objects.filter(
                **{f"{attn}__in": dup_pks}
            ).count()
        elif not rel.related_model._meta.auto_created:
            report.fk_referrers[label] = rel.related_model.objects.filter(
                **{f"{rel.field.name}__in": dup_pks}
            ).count()

    for m2m in Compound._meta.many_to_many:
        through = m2m.remote_field.through
        report.would_fold[f"Compound.{m2m.name}"] = through.objects.filter(
            **{f"{_through_compound_fk_attname(through)}__in": dup_pks}
        ).count()

    _collect_field_conflicts(report)

    return report


def _repoint(dup, keeper):
    """Repoint every referrer of ``dup`` to ``keeper`` (then ``dup`` can be
    deleted).

    NOT ENABLED YET: the mutating statements are commented out until the collapse
    rules are agreed. The reads are left live so the shape is clear; uncomment the
    ``.update()`` / ``.add()`` / ``.delete()`` lines to actually collapse.
    """
    raise NotImplementedError("collapse is disabled pending strategy")

    # model = type(dup)
    #
    # # 1. Concrete reverse FKs, incl. explicit through models (ExperimentCompound,
    # #    CompoundIdentifier). Bulk-repoint; on a unique-constraint clash with a
    # #    row keeper already has, repoint the rest and drop the would-be dupes.
    # for rel in model._meta.related_objects:
    #     if rel.many_to_many or rel.related_model._meta.auto_created:
    #         continue
    #     rel_model, attn = rel.related_model, rel.field.attname
    #     qs = rel_model.objects.filter(**{attn: dup.pk})
    #     try:
    #         with transaction.atomic():
    #             qs.update(**{attn: keeper.pk})
    #     except IntegrityError:
    #         for pk in list(qs.values_list("pk", flat=True)):
    #             try:
    #                 with transaction.atomic():
    #                     rel_model.objects.filter(pk=pk).update(**{attn: keeper.pk})
    #             except IntegrityError:
    #                 rel_model.objects.filter(pk=pk).delete()
    #
    # # 2. Auto-through M2M where Compound is the target (DesignSet.compounds).
    # for rel in model._meta.related_objects:
    #     if not rel.many_to_many or not rel.through._meta.auto_created:
    #         continue
    #     through = rel.through
    #     attn = _through_compound_fk_attname(through)
    #     for pk in list(through.objects.filter(**{attn: dup.pk}).values_list("pk", flat=True)):
    #         try:
    #             with transaction.atomic():
    #                 through.objects.filter(pk=pk).update(**{attn: keeper.pk})
    #         except IntegrityError:
    #             through.objects.filter(pk=pk).delete()
    #
    # # 3. Fold the duplicate's own forward M2M links (inspirations) into keeper.
    # for m2m in model._meta.many_to_many:
    #     getattr(keeper, m2m.name).add(*getattr(dup, m2m.name).all())
