# Generated by Django 3.2.25 on 2024-06-14 10:16

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0057_auto_20240612_1348'),
    ]

    operations = [
        migrations.CreateModel(
            name='ComputedSetComputedMolecule',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('computed_molecule', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.computedmolecule')),
                ('computed_set', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.computedset')),
            ],
        ),
        migrations.AddConstraint(
            model_name='computedsetcomputedmolecule',
            constraint=models.UniqueConstraint(fields=('computed_set', 'computed_molecule'), name='unique_computedsetcomputedmolecule'),
        ),
        migrations.RemoveField(
            model_name='computedmolecule',
            name='computed_set',
        ),
        migrations.AddField(
            model_name='computedset',
            name='computed_molecules',
            field=models.ManyToManyField(related_name='computed_set', through='viewer.ComputedSetComputedMolecule', to='viewer.ComputedMolecule'),
        ),
    ]
