# Generated by Django 3.1.14 on 2022-10-26 09:04

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0027_sessionproject_project'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='protein',
            unique_together={('code', 'target_id', 'prot_type')},
        ),
    ]
