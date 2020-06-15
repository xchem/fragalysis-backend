# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2020-06-12 16:22
from __future__ import unicode_literals

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone
import uuid


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='ActivityPoint',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('source', models.CharField(db_index=True, max_length=50, null=True)),
                ('activity', models.FloatField(db_index=True)),
                ('units', models.CharField(max_length=50)),
                ('confidence', models.IntegerField(db_index=True, null=True)),
                ('internal_id', models.CharField(max_length=150, null=True)),
                ('operator', models.CharField(default=b'NA', max_length=5)),
            ],
            options={
                'permissions': (('view_activitypoint', 'View activitypoint'),),
            },
        ),
        migrations.CreateModel(
            name='Compound',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('inchi', models.CharField(db_index=True, max_length=255, unique=True)),
                ('smiles', models.CharField(db_index=True, max_length=255)),
                ('mol_log_p', models.FloatField()),
                ('mol_wt', models.FloatField()),
                ('tpsa', models.FloatField()),
                ('heavy_atom_count', models.IntegerField()),
                ('heavy_atom_mol_wt', models.FloatField()),
                ('nhoh_count', models.IntegerField()),
                ('no_count', models.IntegerField()),
                ('num_h_acceptors', models.IntegerField()),
                ('num_h_donors', models.IntegerField()),
                ('num_het_atoms', models.IntegerField()),
                ('num_rot_bonds', models.IntegerField()),
                ('num_val_electrons', models.IntegerField()),
                ('ring_count', models.IntegerField()),
            ],
            options={
                'permissions': (('view_compound', 'View compound'),),
            },
        ),
        migrations.CreateModel(
            name='CompoundSet',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50, unique=True)),
                ('submitted_sdf', models.FileField(max_length=255, upload_to=b'compound_sets/')),
                ('spec_version', models.FloatField()),
                ('method_url', models.TextField(max_length=1000, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='CompoundSetSubmitter',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('email', models.CharField(max_length=100)),
                ('institution', models.CharField(max_length=50)),
                ('generation_date', models.DateField()),
                ('method', models.CharField(max_length=50)),
                ('unique_name', models.CharField(max_length=101)),
            ],
        ),
        migrations.CreateModel(
            name='ComputedCompound',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('sdf_info', models.TextField()),
                ('name', models.CharField(max_length=50)),
                ('smiles', models.CharField(max_length=255)),
                ('original_smiles', models.CharField(max_length=255)),
                ('pdb_info', models.FileField(max_length=255, upload_to=b'pdbs/')),
                ('compound_set', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.CompoundSet')),
            ],
        ),
        migrations.CreateModel(
            name='CSetKeys',
            fields=[
                ('user', models.TextField(max_length=50)),
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
            ],
        ),
        migrations.CreateModel(
            name='Molecule',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('smiles', models.CharField(db_index=True, max_length=255, null=True)),
                ('lig_id', models.CharField(max_length=5, null=True)),
                ('chain_id', models.CharField(max_length=1, null=True)),
                ('mol_type', models.CharField(choices=[(b'PR', b'Proasis molecule'), (b'SD', b'Sdf molecule'), (b'HA', b'Hydrogens added '), (b'HC', b'Mol2 format with Hydrogens and AM1 BCC')], default=b'PR', max_length=2)),
                ('sdf_info', models.TextField(null=True)),
                ('rscc', models.FloatField(null=True)),
                ('occupancy', models.FloatField(null=True)),
                ('x_com', models.FloatField(null=True)),
                ('y_com', models.FloatField(null=True)),
                ('z_com', models.FloatField(null=True)),
                ('rmsd', models.FloatField(null=True)),
                ('cmpd_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Compound')),
            ],
            options={
                'permissions': (('view_molecule', 'View molecule'),),
            },
        ),
        migrations.CreateModel(
            name='NumericalScoreValues',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.FloatField()),
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedCompound')),
            ],
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=200, unique=True)),
                ('init_date', models.DateTimeField(auto_now_add=True)),
                ('user_id', models.ManyToManyField(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'permissions': (('view_project', 'View project'),),
            },
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('code', models.CharField(db_index=True, max_length=50)),
                ('apo_holo', models.NullBooleanField()),
                ('prot_type', models.CharField(choices=[(b'AP', b'Apo'), (b'ST', b'Stripped'), (b'TL', b'Tleaped'), (b'CH', b'Chunked'), (b'BO', b'Bound')], default=b'AP', max_length=2)),
                ('pdb_info', models.FileField(max_length=255, null=True, upload_to=b'pdbs/')),
                ('bound_info', models.FileField(max_length=255, null=True, upload_to=b'bound/')),
                ('cif_info', models.FileField(max_length=255, null=True, upload_to=b'cifs/')),
                ('mtz_info', models.FileField(max_length=255, null=True, upload_to=b'mtzs/')),
                ('map_info', models.FileField(max_length=255, null=True, upload_to=b'maps/')),
                ('aligned', models.NullBooleanField()),
                ('has_eds', models.NullBooleanField()),
                ('aligned_to', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.Protein')),
            ],
            options={
                'permissions': (('view_protein', 'View protein'),),
            },
        ),
        migrations.CreateModel(
            name='ScoreDescription',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField()),
                ('compound_set', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.CompoundSet')),
            ],
        ),
        migrations.CreateModel(
            name='SessionProject',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=200)),
                ('init_date', models.DateTimeField(default=django.utils.timezone.now)),
                ('description', models.CharField(default=b'', max_length=255)),
                ('tags', models.TextField(default=b'[]')),
                ('author', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'db_table': 'viewer_sessionproject',
                'permissions': (('view_project', 'View project'),),
            },
        ),
        migrations.CreateModel(
            name='Snapshot',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False)),
                ('type', models.CharField(choices=[(b'INIT', b'INIT'), (b'AUTO', b'AUTO'), (b'MANUAL', b'MANUAL')], default=b'INIT', max_length=8)),
                ('title', models.CharField(max_length=255)),
                ('description', models.CharField(default=b'', max_length=255)),
                ('created', models.DateTimeField(default=django.utils.timezone.now)),
                ('data', models.TextField()),
                ('author', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
                ('parent', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='children', to='viewer.Snapshot')),
                ('session_project', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.SessionProject')),
            ],
            options={
                'db_table': 'viewer_snapshot',
                'managed': True,
                'permissions': (('view_project', 'View project'),),
            },
        ),
        migrations.CreateModel(
            name='Target',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=200, unique=True)),
                ('init_date', models.DateTimeField(auto_now_add=True)),
                ('uniprot_id', models.CharField(max_length=100, null=True)),
                ('project_id', models.ManyToManyField(to='viewer.Project')),
            ],
            options={
                'permissions': (('view_target', 'View target'),),
            },
        ),
        migrations.CreateModel(
            name='TextScoreValues',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.TextField(max_length=500)),
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedCompound')),
                ('score', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ScoreDescription')),
            ],
        ),
        migrations.AddField(
            model_name='sessionproject',
            name='target',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Target'),
        ),
        migrations.AddField(
            model_name='protein',
            name='target_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Target'),
        ),
        migrations.AddField(
            model_name='numericalscorevalues',
            name='score',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ScoreDescription'),
        ),
        migrations.AddField(
            model_name='molecule',
            name='prot_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Protein'),
        ),
        migrations.AddField(
            model_name='computedcompound',
            name='inspiration_frags',
            field=models.ManyToManyField(to='viewer.Molecule'),
        ),
        migrations.AlterUniqueTogether(
            name='compoundsetsubmitter',
            unique_together=set([('name', 'method')]),
        ),
        migrations.AddField(
            model_name='compoundset',
            name='submitter',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.CompoundSetSubmitter'),
        ),
        migrations.AddField(
            model_name='compoundset',
            name='target',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Target'),
        ),
        migrations.AddField(
            model_name='compound',
            name='project_id',
            field=models.ManyToManyField(to='viewer.Project'),
        ),
        migrations.AddField(
            model_name='activitypoint',
            name='cmpd_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Compound'),
        ),
        migrations.AddField(
            model_name='activitypoint',
            name='target_id',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Target'),
        ),
        migrations.AlterUniqueTogether(
            name='protein',
            unique_together=set([('code', 'prot_type')]),
        ),
        migrations.AlterUniqueTogether(
            name='molecule',
            unique_together=set([('prot_id', 'cmpd_id', 'mol_type')]),
        ),
        migrations.AlterUniqueTogether(
            name='activitypoint',
            unique_together=set([('target_id', 'activity', 'cmpd_id', 'units')]),
        ),
    ]
