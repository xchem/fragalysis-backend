"""Investigate / collapse duplicate Compound rows (same project + inchi_key).

Now that ``Compound.project`` is a scalar FK, duplicates are simply the rows
sharing a ``(project, inchi_key)``. Each such group should collapse to exactly one
canonical row, with *every* object that points at a duplicate - through any FK or
M2M - repointed to the survivor.

DRY RUN: :func:`collapse_duplicate_compounds` does NOT modify the database. It
returns a :class:`DedupReport` describing what would happen (which keeper per
group, and how many referring objects each relation would move), so the exact
collapse rules can be decided first. The actual repoint/delete implementation
lives in :func:`_repoint` with its mutating calls commented out, ready to enable.

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
    mergeable_groups: int = 0
    conflict_groups: int = 0
    would_remove: int = 0
    # relation label -> number of referring objects that would move to the keeper
    fk_referrers: dict = field(default_factory=dict)
    m2m_referrers: dict = field(default_factory=dict)
    # forward m2m on Compound whose links would be folded into the keeper
    would_fold: dict = field(default_factory=dict)
    # (project_id, inchi_key, [codes]) for groups left untouched
    conflicts: list = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            "",
            "=== Compound dedup (DRY RUN - nothing written) ===",
            f"  compounds:                    {self.total_compounds}",
            f"  (project, inchi_key) groups:  {self.groups}",
            f"  duplicate groups:             {self.duplicate_groups}",
            f"    mergeable:                  {self.mergeable_groups}"
            f"  -> would remove {self.would_remove} rows",
            f"    conflicts (>1 code):        {self.conflict_groups}  (untouched)",
            "  referrers that would repoint to the keeper (FK):",
        ]
        lines += [f"      {k}: {v}" for k, v in sorted(self.fk_referrers.items())]
        lines.append("  referrers via M2M:")
        lines += [f"      {k}: {v}" for k, v in sorted(self.m2m_referrers.items())]
        lines.append("  keeper absorbs the duplicates' own M2M links:")
        lines += [f"      {k}: {v}" for k, v in sorted(self.would_fold.items())]
        return "\n".join(lines)


def _through_compound_fk_attname(through) -> str:
    """The attname of the FK on ``through`` that points at Compound."""
    for f in through._meta.get_fields():
        if getattr(f, "related_model", None) is Compound and not f.many_to_many:
            return f.attname
    raise LookupError(f"no Compound FK on {through._meta.label}")


def _pick_keeper(members, codes):
    """The surviving row: the coded one (lowest pk) if a single code exists, else
    the lowest-pk row."""
    if codes:
        return min((c for c in members if c.compound_code), key=lambda c: c.pk)
    return min(members, key=lambda c: c.pk)


def collapse_duplicate_compounds() -> DedupReport:
    report = DedupReport()
    report.total_compounds = Compound.objects.count()

    groups: dict = defaultdict(list)
    for cmpd in (
        Compound.objects.exclude(inchi_key="")
        .exclude(inchi_key__isnull=True)
        .only("pk", "project_id", "inchi_key", "compound_code")
    ):
        groups[(cmpd.project_id, cmpd.inchi_key)].append(cmpd)
    report.groups = len(groups)

    dup_pks: list = []
    for members in groups.values():
        if len(members) < 2:
            continue
        report.duplicate_groups += 1

        codes = sorted({c.compound_code for c in members if c.compound_code})
        if len(codes) >= 2:
            report.conflict_groups += 1
            report.conflicts.append(
                (members[0].project_id, members[0].inchi_key, codes)
            )
            continue

        report.mergeable_groups += 1
        keeper = _pick_keeper(members, codes)
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
