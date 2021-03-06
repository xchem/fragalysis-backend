# -*- coding: utf-8 -*-
# Generated by Django 1.11.29 on 2020-07-07 22:23


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
                ('long_inchi', models.TextField(blank=True, max_length=1000, null=True)),
                ('smiles', models.CharField(db_index=True, max_length=255)),
                ('current_identifier', models.CharField(blank=True, db_index=True, max_length=255, null=True)),
                ('all_identifiers', models.TextField(blank=True, null=True)),
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
                ('description', models.TextField(blank=True, null=True)),
                ('comments', models.TextField(blank=True, null=True)),
            ],
            options={
                'permissions': (('view_compound', 'View compound'),),
            },
        ),
        migrations.CreateModel(
            name='ComputedMolecule',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('sdf_info', models.TextField()),
                ('name', models.CharField(max_length=50)),
                ('smiles', models.CharField(max_length=255)),
                ('pdb_info', models.FileField(max_length=255, upload_to=b'pdbs/')),
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='ComputedSet',
            fields=[
                ('name', models.CharField(max_length=50, primary_key=True, serialize=False, unique=True)),
                ('submitted_sdf', models.FileField(max_length=255, upload_to=b'compound_sets/')),
                ('spec_version', models.FloatField()),
                ('method_url', models.TextField(max_length=1000, null=True)),
                ('unique_name', models.CharField(max_length=101)),
            ],
        ),
        migrations.CreateModel(
            name='ComputedSetSubmitter',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
                ('email', models.CharField(max_length=100)),
                ('institution', models.CharField(max_length=50)),
                ('generation_date', models.DateField()),
                ('method', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='CSetKeys',
            fields=[
                ('user', models.CharField(default=b'User', editable=False, max_length=50)),
                ('uuid', models.UUIDField(default=uuid.uuid4, editable=False, primary_key=True, serialize=False)),
            ],
        ),
        migrations.CreateModel(
            name='DesignSet',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('set_name', models.CharField(max_length=50)),
                ('set_type', models.CharField(choices=[(b'library', b'library'), (b'follow-up', b'follow-up'), (b'user-submitted', b'user-submitted'), (b'enumerated', b'enumerated')], default=b'user-submitted', max_length=100)),
                ('set_description', models.TextField(blank=True, max_length=1000, null=True)),
                ('compounds', models.ManyToManyField(to='viewer.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='File',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file', models.FileField(upload_to=b'')),
            ],
        ),
        migrations.CreateModel(
            name='HistoricalMolecule',
            fields=[
                ('id', models.IntegerField(auto_created=True, blank=True, db_index=True, verbose_name='ID')),
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
                ('history_id', models.AutoField(primary_key=True, serialize=False)),
                ('history_date', models.DateTimeField()),
                ('history_change_reason', models.CharField(max_length=100, null=True)),
                ('history_type', models.CharField(choices=[('+', 'Created'), ('~', 'Changed'), ('-', 'Deleted')], max_length=1)),
                ('cmpd_id', models.ForeignKey(blank=True, db_constraint=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='+', to='viewer.Compound')),
                ('history_user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='+', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ('-history_date', '-history_id'),
                'get_latest_by': 'history_date',
                'verbose_name': 'historical molecule',
            },
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
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedMolecule')),
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
                ('computed_set', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedSet')),
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
                ('metadata', models.FileField(max_length=255, null=True, upload_to=b'metadata/')),
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
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedMolecule')),
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
            model_name='historicalmolecule',
            name='prot_id',
            field=models.ForeignKey(blank=True, db_constraint=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='+', to='viewer.Protein'),
        ),
        migrations.AlterUniqueTogether(
            name='computedsetsubmitter',
            unique_together=set([('name', 'method')]),
        ),
        migrations.AddField(
            model_name='computedset',
            name='submitter',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedSetSubmitter'),
        ),
        migrations.AddField(
            model_name='computedset',
            name='target',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.Target'),
        ),
        migrations.AddField(
            model_name='computedmolecule',
            name='computed_inspirations',
            field=models.ManyToManyField(blank=True, null=True, to='viewer.Molecule'),
        ),
        migrations.AddField(
            model_name='computedmolecule',
            name='computed_set',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.ComputedSet'),
        ),
        migrations.AddField(
            model_name='compound',
            name='inspirations',
            field=models.ManyToManyField(blank=True, null=True, to='viewer.Molecule'),
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
