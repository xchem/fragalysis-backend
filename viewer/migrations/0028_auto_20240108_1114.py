# Generated by Django 3.2.23 on 2024-01-08 11:14

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0027_rename_lhs_pdb_in_computedmolecule'),
    ]

    operations = [
        migrations.RenameField(
            model_name='historicalsiteobservation',
            old_name='xmap_fofc_file',
            new_name='diff_file',
        ),
        migrations.RenameField(
            model_name='historicalsiteobservation',
            old_name='xmap_2fofc_file',
            new_name='sigmaa_file',
        ),
        migrations.RenameField(
            model_name='siteobservation',
            old_name='xmap_fofc_file',
            new_name='diff_file',
        ),
        migrations.RenameField(
            model_name='siteobservation',
            old_name='xmap_2fofc_file',
            new_name='sigmaa_file',
        ),
    ]
