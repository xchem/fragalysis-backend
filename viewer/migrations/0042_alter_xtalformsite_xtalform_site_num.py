# Generated by Django 3.2.23 on 2024-02-05 12:02

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0041_xtalform_xtalform_num'),
    ]

    operations = [
        migrations.AlterField(
            model_name='xtalformsite',
            name='xtalform_site_num',
            field=models.TextField(
                help_text='alphabetic xtalform site id (enumerated on creation)',
                null=True,
            ),
        ),
    ]
