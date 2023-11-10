# Generated by Django 3.2.20 on 2023-11-10 11:51

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0017_merge_20231109_1252'),
    ]

    operations = [
        migrations.AddField(
            model_name='historicalsiteobservation',
            name='xmap_fofc_file',
            field=models.TextField(max_length=255, null=True),
        ),
        migrations.AddField(
            model_name='siteobservation',
            name='xmap_fofc_file',
            field=models.FileField(max_length=255, null=True, upload_to='target_loader_data/'),
        ),
    ]
