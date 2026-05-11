"""Convert DownloadLinks.file_url from absolute path to basename.

The absolute path is now reconstructed at access time via
DownloadLinks.get_file_url() using MEDIA_ROOT, the "downloads" subdir and
task_id. This migration:

  1. Drops the unique=True constraint, since two records for different tasks
     can legitimately produce the same filename (e.g. "TARGET.zip"). This
     has to happen first — otherwise the backfill below would hit the
     constraint as soon as it tries to save a duplicate basename.
  2. Backfills existing rows: copies the per-download UUID directory from
     the old absolute path into task_id and replaces
     file_url with the original basename.
"""
import os

from django.db import migrations, models


def _strip_path(apps, schema_editor):
    DownloadLinks = apps.get_model('viewer', 'DownloadLinks')
    for record in DownloadLinks.objects.exclude(file_url__isnull=True).exclude(
        file_url__exact=''
    ):
        old_value = record.file_url
        # Old format was "{MEDIA_ROOT}/downloads/{uuid}/{filename}".
        # Replace old task_ids with the parent directory name (the original uuid) if present.
        record.task_id = os.path.basename(os.path.dirname(old_value)) or record.task_id
        record.file_url = os.path.basename(old_value)
        record.save(update_fields=['file_url', 'task_id'])


def _noop_reverse(apps, schema_editor):
    # Reversing the data change isn't safe — we'd have to guess MEDIA_ROOT.
    # Schema rollback is handled by the AlterField below; data stays as-is.
    pass


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0154_project_alias_alter_project_title'),
    ]

    operations = [
        migrations.AlterField(
            model_name='downloadlinks',
            name='file_url',
            field=models.TextField(db_index=True, null=True),
        ),
        migrations.RunPython(_strip_path, _noop_reverse),
    ]
