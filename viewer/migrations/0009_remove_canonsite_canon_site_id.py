# Generated by Django 3.2.20 on 2023-08-29 20:16

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0008_alter_experiment_xtalform'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='canonsite',
            name='canon_site_id',
        ),
    ]