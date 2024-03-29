# Generated by Django 3.2.23 on 2024-01-11 12:03

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0029_target_display_name'),
    ]

    operations = [
        migrations.AddField(
            model_name='computedmolecule',
            name='reference_code',
            field=models.TextField(
                blank=True,
                help_text='The computed reference SiteObservation (the corresponding ref_pdb value if it has one)',
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name='computedmolecule',
            name='site_observation_code',
            field=models.TextField(
                blank=True,
                help_text='The LHS SiteObservation (the corresponding lhs_pdb value if it has one)',
                null=True,
            ),
        ),
    ]
