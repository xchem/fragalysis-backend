# Generated by Django 3.1 on 2022-01-11 15:51

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('scoring', '0003_auto_20200804_1411'),
        ('viewer', '0018_auto_20211214_1609'),
    ]

    operations = [
        migrations.AlterField(
            model_name='moleculetag',
            name='mol_group',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='scoring.molgroup'),
        ),
    ]