"""Step 2/3 of moving Compound.project_id (M2M) to Compound.project (FK).

Backfills the FK from the legacy M2M, choosing a single project per compound:

* A compound's SiteObservations are authoritative - their target's project is the
  real one (SiteObservation -> experiment -> experiment_upload -> target ->
  project). If the observations name exactly one project, that wins.
* The legacy ``project_id`` M2M is the fallback for compounds that have no
  observations (e.g. computed-set compounds).
* Ambiguities - a compound observed across more than one project, or a
  multi-project M2M with no observations - are resolved deterministically to the
  lowest project id and reported.
* Compounds with neither an observation nor an M2M link can't be assigned; they
  are reported and left null, and will block 0168 (NOT NULL) until resolved.

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


def backfill(apps, schema_editor):
    Compound = apps.get_model("viewer", "Compound")
    SiteObservation = apps.get_model("viewer", "SiteObservation")
    Project = apps.get_model("viewer", "Project")

    # M2M links, read via the through table (unambiguous).
    through = Compound._meta.get_field("project_id").remote_field.through
    m2m = defaultdict(set)
    for compound_id, project_id in through.objects.values_list(
        "compound_id", "project_id"
    ):
        m2m[compound_id].add(project_id)

    # Authoritative project(s) from each compound's site observations.
    so = defaultdict(set)
    for compound_id, project_id in (
        SiteObservation.objects.filter(cmpd__isnull=False)
        .values_list("cmpd_id", "experiment__experiment_upload__target__project")
        .distinct()
    ):
        if project_id is not None:
            so[compound_id].add(project_id)

    projects = {p.pk: p for p in Project.objects.all()}

    by_project: dict = defaultdict(list)  # project_id -> [compound_id, ...]
    unresolved: list = []
    stats = {
        "so_single": 0,
        "so_disagree": 0,
        "so_multi": 0,
        "m2m_single": 0,
        "m2m_multi": 0,
    }

    for compound_id in Compound.objects.values_list("pk", flat=True):
        obs_projects = so.get(compound_id, set())
        m2m_projects = m2m.get(compound_id, set())

        if len(obs_projects) == 1:
            chosen = next(iter(obs_projects))
            stats["so_single"] += 1
            if m2m_projects and m2m_projects != obs_projects:
                stats["so_disagree"] += 1
        elif len(obs_projects) > 1:
            chosen = min(obs_projects)
            stats["so_multi"] += 1
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
        f"  from site observations, 1 project:   {stats['so_single']}"
        f"  (disagreed with M2M: {stats['so_disagree']})",
        f"  observed across >1 project (min id):  {stats['so_multi']}",
        f"  no observations, single M2M project:  {stats['m2m_single']}",
        f"  no observations, multi M2M (min id):  {stats['m2m_multi']}",
        f"  UNRESOLVED (no observation, no M2M):  {len(unresolved)}",
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
