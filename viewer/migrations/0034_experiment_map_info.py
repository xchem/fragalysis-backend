# Generated by Django 3.2.23 on 2024-01-26 11:16

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0033_alter_siteobservation_cmpd'),
    ]

    operations = [
        migrations.AddField(
            model_name='experiment',
            name='map_info',
            field=models.FileField(
                max_length=255, null=True, upload_to='target_loader_data/'
            ),
        ),
    ]
