# Generated by Django 3.2.20 on 2023-08-18 15:07

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    dependencies = [
        ('scoring', '0001_initial'),
        ('viewer', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='siteobservationtag',
            name='mol_group',
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.SET_NULL,
                to='scoring.siteobservationgroup',
            ),
        ),
    ]
