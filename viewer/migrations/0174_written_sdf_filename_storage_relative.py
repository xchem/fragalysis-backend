"""Make ComputedSet.written_sdf_filename a storage name (issue 1038).

The computed-set upload used to record an absolute path under MEDIA_ROOT
(i.e. "/code/media/computed_set_data/x.sdf"). With S3 media storage there is
no MEDIA_ROOT, so the value is now a storage name ("computed_set_data/x.sdf"),
like a FileField's name. This changes the existing rows, in the model's table
and in its history table (the computed-set zip downloads read the history).

The reverse migration puts MEDIA_ROOT back in front of every storage name.
"""
from django.conf import settings
from django.db import migrations
from django.db.models import F, Value
from django.db.models.functions import Concat, Substr

_MODELS = ("ComputedSet", "HistoricalComputedSet")


def _media_root() -> str:
    return str(settings.MEDIA_ROOT).rstrip("/") + "/"


def to_storage_names(apps, schema_editor):
    root = _media_root()
    for model_name in _MODELS:
        model = apps.get_model("viewer", model_name)
        model.objects.filter(written_sdf_filename__startswith=root).update(
            written_sdf_filename=Substr("written_sdf_filename", len(root) + 1)
        )


def to_absolute_paths(apps, schema_editor):
    root = _media_root()
    for model_name in _MODELS:
        model = apps.get_model("viewer", model_name)
        model.objects.exclude(written_sdf_filename__isnull=True).exclude(
            written_sdf_filename=""
        ).exclude(written_sdf_filename__startswith="/").update(
            written_sdf_filename=Concat(Value(root), F("written_sdf_filename"))
        )


class Migration(migrations.Migration):
    dependencies = [
        ("viewer", "0173_target_blog"),
    ]

    operations = [
        migrations.RunPython(to_storage_names, to_absolute_paths),
    ]
