"""Step 3/3 of moving Compound.project_id (M2M) to Compound.project (FK).

Drops the old ``project_id`` M2M - the point of no return for the M2M links.

0166 backfills ``project`` from this M2M, so ordinarily migrating straight to
head is fine. The one thing to check first: 0166 prints an UNRESOLVED count for
compounds it could not assign a project to. If that is non-zero, resolve those
compounds before applying this (once the M2M is gone they can't be backfilled,
and 0168's NOT NULL will fail on them). To inspect first:

    python manage.py migrate viewer 0166   # then read the report
"""

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("viewer", "0166_backfill_compound_project"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="compound",
            name="project_id",
        ),
    ]
