# Generated by Django 3.2.23 on 2024-01-30 08:09

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0034_experiment_map_info'),
    ]

    operations = [
        migrations.AlterField(
            model_name='experiment',
            name='event_map_info',
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.FileField(max_length=255, upload_to=''),
                null=True,
                size=None,
            ),
        ),
    ]