# Generated by Django 3.2.20 on 2023-10-23 10:49

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0014_auto_20231023_1048'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='siteobservationtag',
            name='site_observations',
        ),
        migrations.AddField(
            model_name='siteobservationtag',
            name='site_observations',
            field=models.ManyToManyField(through='viewer.SiteObvsSiteObservationTag', to='viewer.SiteObservation'),
        ),
    ]
