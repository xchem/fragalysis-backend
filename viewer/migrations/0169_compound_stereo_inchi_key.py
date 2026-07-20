"""Recompute Compound.inchi_key to be stereo-inclusive.

Historically the loader computed ``inchi_key`` from a *flattened* molecule
(``Chem.RemoveStereochemistry`` before ``MolToInchiKey``), so enantiomers and
diastereomers collapsed to the same key. That is wrong: they are distinct
compounds. The code that produced flattened keys has been removed; this data
migration brings existing rows in line by recomputing the key from each
compound's ``smiles`` *without* flattening.

No data is lost:

* Target-loaded compounds kept full stereochemistry in ``smiles``, so their key
  is upgraded flat -> stereo-inclusive.
* Compounds whose ``smiles`` was itself stored flat (older computed-set uploads)
  simply recompute to the same flat key - never worse than before.
* Rows with empty/unparseable ``smiles`` are left untouched and reported.

The reverse operation re-flattens, restoring the previous keys exactly for the
rows this migration changed.
"""

from django.db import migrations
from rdkit import Chem


def _log(schema_editor, message):
    writer = getattr(schema_editor.connection, "_migration_stdout", None)
    if writer is not None:
        writer.write(message + "\n")
    else:
        print(message)


def _recompute(apps, schema_editor, flatten):
    Compound = apps.get_model("viewer", "Compound")

    changed = []
    skipped = 0
    for pk, smiles, old_key in Compound.objects.values_list(
        "pk", "smiles", "inchi_key"
    ).iterator():
        if not smiles:
            skipped += 1
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            skipped += 1
            continue
        if flatten:
            Chem.RemoveStereochemistry(mol)
        new_key = Chem.inchi.MolToInchiKey(mol)
        if new_key and new_key != old_key:
            changed.append((pk, new_key))

    # Bulk update in chunks to keep memory/SQL sizes reasonable.
    Compound.objects.bulk_update(
        [Compound(pk=pk, inchi_key=key) for pk, key in changed],
        ["inchi_key"],
        batch_size=1000,
    )

    _log(
        schema_editor,
        f"inchi_key recompute ({'flatten' if flatten else 'stereo'}): "
        f"{len(changed)} updated, {skipped} skipped (empty/unparseable smiles)",
    )


def forwards(apps, schema_editor):
    _recompute(apps, schema_editor, flatten=False)


def backwards(apps, schema_editor):
    _recompute(apps, schema_editor, flatten=True)


class Migration(migrations.Migration):

    dependencies = [
        ("viewer", "0168_compound_project_not_null"),
    ]

    operations = [
        migrations.RunPython(forwards, backwards),
    ]
