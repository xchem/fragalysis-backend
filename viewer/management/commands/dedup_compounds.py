"""Report the duplicate-Compound collapse plan (dry run - nothing is written).

Groups compounds by ``(project, inchi_key)``, works out the single keeper per
group, and prints every object that would be repointed to it (via FK or M2M), so
the collapse strategy can be decided before anything is committed. See
``viewer.compound_dedup``.
"""

from django.core.management.base import BaseCommand

from viewer.compound_dedup import collapse_duplicate_compounds


class Command(BaseCommand):
    help = "Dry-run report of the duplicate-compound collapse (writes nothing)."

    def add_arguments(self, parser):
        parser.add_argument(
            "--conflicts",
            action="store_true",
            help="Also list the conflict groups (same structure, >1 code).",
        )

    def handle(self, *args, **kwargs):
        del args
        report = collapse_duplicate_compounds()
        self.stdout.write(report.summary())

        if kwargs["conflicts"] and report.conflicts:
            self.stdout.write("")
            self.stdout.write("  conflict groups (project_id, inchi_key, codes):")
            for project_id, inchi_key, codes in report.conflicts[:100]:
                self.stdout.write(f"    project={project_id} {inchi_key}: {codes}")
            if len(report.conflicts) > 100:
                self.stdout.write(f"    ... and {len(report.conflicts) - 100} more")
