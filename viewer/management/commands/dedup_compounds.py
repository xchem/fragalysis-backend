"""Report the duplicate-Compound collapse plan (dry run - nothing is written).

    ANALYSIS TOOL - NOT PART OF THE APPLICATION.
    Read-only reporting command; it writes nothing. Deduplication was decided to
    require human curation, so this is not wired into any app flow. Safe to
    delete this command and ``viewer/compound_dedup.py`` at any time.

Groups compounds by ``(project, inchi_key, smiles, compound_code)``, works out
the single keeper per group (lowest pk), and prints every object that would be
repointed to it (via FK or M2M), plus content-field conflicts within those
groups, to inform manual curation. See ``viewer.compound_dedup``.
"""

from django.core.management.base import BaseCommand

from viewer.compound_dedup import collapse_duplicate_compounds


class Command(BaseCommand):
    help = "Dry-run report of the duplicate-compound collapse (writes nothing)."

    def add_arguments(self, parser):
        parser.add_argument(
            "--groups",
            action="store_true",
            help="Also list a sample of the duplicate groups that would collapse.",
        )
        parser.add_argument(
            "--conflicts",
            action="store_true",
            help="Also list groups sharing (project, inchi_key, smiles) but with "
            "conflicting values in other content fields.",
        )

    def handle(self, *args, **kwargs):
        del args
        report = collapse_duplicate_compounds()
        self.stdout.write(report.summary())

        if kwargs["groups"] and report.sample_groups:
            self.stdout.write("")
            self.stdout.write(
                "  duplicate groups (project_id, inchi_key, smiles, code, count):"
            )
            for project_id, inchi_key, smiles, code, count in report.sample_groups[
                :100
            ]:
                self.stdout.write(
                    f"    project={project_id} {inchi_key} {smiles!r} "
                    f"code={code!r}: {count} rows"
                )
            if len(report.sample_groups) > 100:
                self.stdout.write(f"    ... and {len(report.sample_groups) - 100} more")

        if kwargs["conflicts"] and report.field_conflicts:
            self.stdout.write("")
            self.stdout.write(
                "  field conflicts within "
                "(project, inchi_key, smiles, compound_code) groups:"
            )
            for (
                project_id,
                inchi_key,
                smiles,
                code,
                conflicts,
            ) in report.field_conflicts[:100]:
                self.stdout.write(
                    f"    project={project_id} {inchi_key} {smiles!r} code={code!r}"
                )
                for fieldname, values in sorted(conflicts.items()):
                    shown = ", ".join(f"{str(v)[:60]!r}" for v in values)
                    self.stdout.write(f"        {fieldname}: {shown}")
            if len(report.field_conflicts) > 100:
                self.stdout.write(
                    f"    ... and {len(report.field_conflicts) - 100} more"
                )
