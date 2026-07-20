import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    """Step 1/3 of moving Compound.project_id (M2M) to Compound.project (FK).

    Adds the new nullable FK alongside the existing M2M. The data migration
    (0166) backfills it from the M2M, and 0167 drops the M2M. Splitting it this
    way keeps the M2M readable while the backfill runs.
    """

    dependencies = [
        ("viewer", "0164_remove_textscorevalues_compound_and_more"),
    ]

    operations = [
        migrations.AddField(
            model_name="compound",
            name="project",
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                to="viewer.project",
            ),
        ),
    ]
