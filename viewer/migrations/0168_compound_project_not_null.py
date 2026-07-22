"""Final step of the Compound.project move: enforce NOT NULL.

The FK was added nullable (0165) so it could be attached to existing rows before
the backfill. Once every compound has a project (0166 backfill), tighten it.

This fails if any compound is still null - i.e. if 0166 reported UNRESOLVED
compounds that weren't dealt with. On a fresh/empty database (e.g. the test DB)
it applies trivially.
"""

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("viewer", "0167_remove_compound_project_id"),
    ]

    operations = [
        migrations.AlterField(
            model_name="compound",
            name="project",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                to="viewer.project",
            ),
        ),
    ]
