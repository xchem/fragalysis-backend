"""Drop the orphaned unique constraint on DownloadLinks.file_url and move its
index into Meta.indexes.

Background - the unique constraint has a see-saw history:

  * 0001_initial   created file_url with unique=True
  * 0155           deliberately removed unique=True (file_url now holds only a
                   basename; two different tasks can legitimately share a
                   filename such as "TARGET.zip" - uniqueness lives in task_id)
  * 0160           was generated from a stale model state (it predates the 0155
                   branch) and its AlterField silently re-added unique=True

Because of the branch/merge around 0163, the migration graph's *topological*
order is 0160 (unique=True) -> 0151 -> 0155 (unique=False) -> here. So Django's
migration STATE already believes file_url is non-unique, while the actual
DATABASE still carries the unique constraint that 0160 applied at runtime. The
two have diverged, which is why a plain AlterField(unique=False) does nothing -
Django compares against its state (already non-unique), sees no change and
emits no SQL, leaving the database constraint in place.

Steps:

  1. SeparateDatabaseAndState realigns the migration STATE to match what the
     database really has (db_index=True, unique=True), emitting NO SQL.
  2. AlterField drops both unique and the field-level index: Django issues a
     real DROP CONSTRAINT (the hashed name is found by introspection, so this
     is portable across databases).
  3. AddIndex re-creates the index the model needs, now declared the modern way
     in Meta.indexes (db_index=True on the field is the old convention).

After this, model, migration state and database all agree, and
``makemigrations`` reports no changes.
"""

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0164_remove_textscorevalues_compound_and_more'),
    ]

    operations = [
        # 1. State-only: tell Django the field IS unique (matching the real DB).
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.AlterField(
                    model_name='downloadlinks',
                    name='file_url',
                    field=models.TextField(db_index=True, null=True, unique=True),
                ),
            ],
            database_operations=[],
        ),
        # 2. Real change: drop unique (and the field-level index). Django finds
        #    the unique constraint by introspection and issues DROP CONSTRAINT.
        migrations.AlterField(
            model_name='downloadlinks',
            name='file_url',
            field=models.TextField(null=True),
        ),
        # 3. Re-add the index the modern way, via Meta.indexes.
        migrations.AddIndex(
            model_name='downloadlinks',
            index=models.Index(
                fields=['file_url'], name='downloadlinks_file_url_idx'
            ),
        ),
    ]
