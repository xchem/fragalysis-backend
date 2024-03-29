# Generated by Django 3.2.23 on 2024-02-02 11:35

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0038_auto_20240202_1030'),
    ]

    operations = [
        migrations.AddField(
            model_name='canonsite',
            name='canon_site_num',
            field=models.IntegerField(
                help_text='numeric canon site id (enumerated on creation)', null=True
            ),
        ),
    ]
