# Generated by Django 3.2.23 on 2024-01-30 08:12

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0036_remove_experiment_map_info'),
    ]

    operations = [
        migrations.RenameField(
            model_name='experiment',
            old_name='event_map_info',
            new_name='map_info',
        ),
    ]