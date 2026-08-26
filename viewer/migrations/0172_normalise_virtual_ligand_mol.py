from pathlib import Path, PurePosixPath

from django.conf import settings
from django.db import migrations


def normalise_virtual_ligand_mol(apps, schema_editor):
    """Give every ``virtual_ligand_mol`` value the same shape.

    ``SiteObservation.virtual_ligand_mol`` holds a path relative to
    MEDIA_ROOT - that is the single convention, and what
    ``viewer.cset_upload`` has written since commit 726624a4.

    Migration 0149 wrote a bare basename instead, with the file itself in
    ``computed_set_data/``. Two conventions in one column meant no caller
    could resolve it without first guessing which one it was looking at,
    and ``viewer.download_structures`` guessed wrong for everything
    cset_upload had written - see issue #1025.

    Rewriting those basenames to ``computed_set_data/<name>`` points them at
    the files they already refer to. No files are moved; only the column is
    corrected to say where they are.
    """
    SiteObservation = apps.get_model("viewer", "SiteObservation")

    computed_set_dir = settings.COMPUTED_SET_MEDIA_DIRECTORY
    updated = []

    # Values that already carry a directory are in the target convention
    # (whether computed_set_data/... or target_loader_data/...) and are left
    # alone, which also makes this migration idempotent.
    for pk, value in (
        SiteObservation.objects.exclude(virtual_ligand_mol__isnull=True)
        .exclude(virtual_ligand_mol="")
        .values_list("pk", "virtual_ligand_mol")
    ):
        name = str(value).strip()
        if not name or PurePosixPath(name).parent != PurePosixPath("."):
            continue

        updated.append(
            SiteObservation(
                pk=pk,
                virtual_ligand_mol=str(Path(computed_set_dir).joinpath(name)),
            )
        )

    if updated:
        SiteObservation.objects.bulk_update(
            updated, ["virtual_ligand_mol"], batch_size=1000
        )


def reverse(apps, schema_editor):
    """Deliberately a no-op.

    Reversing would mean reintroducing the bare-basename convention this
    migration exists to remove, and the forward direction is not lossy: a
    ``computed_set_data/<name>`` value names exactly the file its basename
    used to name.
    """


class Migration(migrations.Migration):
    dependencies = [
        ("viewer", "0171_resolve_duplicate_compounds"),
    ]

    operations = [
        migrations.RunPython(normalise_virtual_ligand_mol, reverse),
    ]
