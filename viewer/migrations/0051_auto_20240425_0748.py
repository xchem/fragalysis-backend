# Generated by Django 3.2.24 on 2024-04-25 07:48

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0050_auto_20240412_0930'),
    ]

    operations = [
        migrations.AddField(
            model_name='historicalsiteobservation',
            name='version',
            field=models.PositiveSmallIntegerField(default=1),
        ),
        migrations.AddField(
            model_name='siteobservation',
            name='version',
            field=models.PositiveSmallIntegerField(default=1),
        ),
    ]
