# Generated by Django 3.2.20 on 2023-08-29 07:10

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0005_auto_20230829_0659'),
    ]

    operations = [
        migrations.CreateModel(
            name='ExperimentCompound',
            fields=[
                (
                    'id',
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name='ID',
                    ),
                ),
                (
                    'compound',
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        to='viewer.compound',
                    ),
                ),
                (
                    'experiment',
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        to='viewer.experiment',
                    ),
                ),
            ],
        ),
        migrations.AddField(
            model_name='experiment',
            name='compounds',
            field=models.ManyToManyField(
                through='viewer.ExperimentCompound', to='viewer.Compound'
            ),
        ),
        migrations.AddConstraint(
            model_name='experimentcompound',
            constraint=models.UniqueConstraint(
                fields=('experiment', 'compound'), name='unique_experimentcompound'
            ),
        ),
    ]
