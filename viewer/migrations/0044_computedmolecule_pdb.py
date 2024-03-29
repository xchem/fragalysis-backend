# Generated by Django 3.2.23 on 2024-02-20 15:09

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0043_experiment_prefix_tooltip'),
    ]

    operations = [
        migrations.AddField(
            model_name='computedmolecule',
            name='pdb',
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name='pdb',
                to='viewer.siteobservation',
            ),
        ),
    ]
