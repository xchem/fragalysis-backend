# Generated by Django 3.2.20 on 2023-08-30 09:43

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0009_remove_canonsite_canon_site_id'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='xtalform',
            name='xtalform_id',
        ),
    ]
