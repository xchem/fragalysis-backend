# Generated by Django 3.2.24 on 2024-04-10 15:30

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0050_auto_20240410_1125'),
    ]

    operations = [
        migrations.AlterField(
            model_name='siteobservation',
            name='pose',
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name='site_observations',
                to='viewer.pose',
            ),
        ),
    ]
