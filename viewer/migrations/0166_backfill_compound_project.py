"""Step 2/3 of moving Compound.project_id (M2M) to Compound.project (FK).

Backfills the scalar FK by choosing a single project per compound. A compound is
tied to a project through several structural paths, all rooted in the upload that
created the rows; these are authoritative and are unioned together:

* SiteObservation   -> experiment -> experiment_upload -> project
* ExperimentCompound -> experiment -> experiment_upload -> project
* Result (computed) -> computed_set -> target -> project

NB on the computed-set path: for data that predates this migration it yields
nothing, because the old ``ComputedMolecule`` link (compound -> computed set)
is dropped by 0164 *before* this backfill runs, and that deletion is not
migrated into ``Result.computed_set``. Those computed-set compounds are
therefore covered by the M2M fallback below (verified against real data: the
M2M project matches the computed set's project for every such compound). The
``Result.computed_set`` derivation is kept because it is correct for any data
where that link *is* populated (e.g. newer uploads / re-runs).

The legacy ``project_id`` M2M is only a *fallback*, used when a compound has no
structural link at all (e.g. compounds that exist solely as a project-scoped
row). Precedence:

  1. exactly one structural project            -> use it
  2. several structural projects (cross-project) -> lowest id, reported
  3. no structural link, one M2M project        -> use it
  4. no structural link, several M2M projects   -> lowest id, reported
  5. nothing at all                             -> unresolved, reported

Ambiguities (2/4) are resolved deterministically to the lowest project id so the
migration is repeatable. Unresolved compounds (5) are left null and will block
0168 (NOT NULL) until dealt with.

Why the structural paths and not just the M2M: the M2M is a derived cache and
may be incomplete. Deriving from the upload relations directly means a compound
is never left unresolved (or assigned a stale project) when a hard structural
link exists. Accuracy is the goal here, not speed - everything is pulled into
memory and compared explicitly.

Implementation note: while the M2M ``project_id`` and the new FK ``project``
(whose attname is also ``project_id``) coexist, that name is ambiguous. So we
read the M2M through its *through table* and write the FK with ``.update()``,
neither of which touches the ambiguous instance/field name.
"""

from collections import defaultdict

from django.db import migrations


def _log(schema_editor, message):
    writer = getattr(schema_editor.connection, "_migration_stdout", None)
    if writer is not None:
        writer.write(message + "\n")
    else:
        print(message)


def _group(pairs):
    """(entity_id, project_id) pairs -> {entity_id: {project_id, ...}}, skipping
    rows where the project came out NULL (e.g. computed site observations with no
    experiment, or a computed set with no target)."""
    out: dict = defaultdict(set)
    for entity_id, project_id in pairs:
        if project_id is not None:
            out[entity_id].add(project_id)
    return out


def backfill(apps, schema_editor):
    Compound = apps.get_model("viewer", "Compound")
    SiteObservation = apps.get_model("viewer", "SiteObservation")
    ExperimentCompound = apps.get_model("viewer", "ExperimentCompound")
    Result = apps.get_model("viewer", "Result")
    Project = apps.get_model("viewer", "Project")

    # --- Structural sources (authoritative), all rooted in the upload ---------

    # 1. SiteObservation -> experiment -> experiment_upload -> project
    so = _group(
        SiteObservation.objects.filter(cmpd__isnull=False)
        .values_list("cmpd_id", "experiment__experiment_upload__project")
        .distinct()
    )

    # 2. ExperimentCompound -> experiment -> experiment_upload -> project
    ec = _group(
        ExperimentCompound.objects.values_list(
            "compound_id", "experiment__experiment_upload__project"
        ).distinct()
    )

    # 3. Result (computed) -> computed_set -> target -> project
    res = _group(
        Result.objects.filter(
            compound__isnull=False, computed_set__isnull=False
        )
        .values_list("compound_id", "computed_set__target__project")
        .distinct()
    )

    # --- Legacy M2M (fallback only), read via the through table (unambiguous) --
    through = Compound._meta.get_field("project_id").remote_field.through
    m2m = _group(
        through.objects.values_list("compound_id", "project_id")
    )

    projects = {p.pk: p for p in Project.objects.all()}

    by_project: dict = defaultdict(list)  # project_id -> [compound_id, ...]
    unresolved: list = []
    stats = {
        "structural_single": 0,
        "structural_multi": 0,
        "structural_disagreed_m2m": 0,
        "m2m_single": 0,
        "m2m_multi": 0,
    }

    for compound_id in Compound.objects.values_list("pk", flat=True):
        structural = (
            so.get(compound_id, set())
            | ec.get(compound_id, set())
            | res.get(compound_id, set())
        )
        m2m_projects = m2m.get(compound_id, set())

        if len(structural) >= 1:
            chosen = min(structural)
            if len(structural) == 1:
                stats["structural_single"] += 1
            else:
                stats["structural_multi"] += 1
            if m2m_projects and m2m_projects != structural:
                stats["structural_disagreed_m2m"] += 1
        elif len(m2m_projects) == 1:
            chosen = next(iter(m2m_projects))
            stats["m2m_single"] += 1
        elif len(m2m_projects) > 1:
            chosen = min(m2m_projects)
            stats["m2m_multi"] += 1
        else:
            unresolved.append(compound_id)
            continue

        by_project[chosen].append(compound_id)

    for project_id, compound_ids in by_project.items():
        Compound.objects.filter(pk__in=compound_ids).update(
            project=projects[project_id]
        )

    lines = [
        "",
        "=== Compound.project backfill ===",
        f"  structural, 1 project:                {stats['structural_single']}"
        f"  (disagreed with M2M: {stats['structural_disagreed_m2m']})",
        f"  structural, >1 project (min id):      {stats['structural_multi']}",
        f"  no structural link, single M2M:       {stats['m2m_single']}",
        f"  no structural link, multi M2M (min):  {stats['m2m_multi']}",
        f"  UNRESOLVED (no structural, no M2M):   {len(unresolved)}",
    ]
    if unresolved:
        lines.append(f"    unresolved ids (first 50): {sorted(unresolved)[:50]}")
        lines.append("    ^ resolve these before applying 0168 (NOT NULL)")
    _log(schema_editor, "\n".join(lines))


class Migration(migrations.Migration):

    dependencies = [
        ("viewer", "0165_compound_project"),
    ]

    operations = [
        # Reverse is a no-op: we don't restore the M2M links.
        migrations.RunPython(backfill, migrations.RunPython.noop),
    ]
