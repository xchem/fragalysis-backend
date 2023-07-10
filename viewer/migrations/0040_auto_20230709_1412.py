# Generated by Django 3.1.14 on 2023-07-09 14:12

from django.conf import settings
import django.core.serializers.json
from django.db import migrations, models
import django.db.models.deletion
import simple_history.models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('viewer', '0039_auto_20230709_1308'),
    ]

    operations = [
        migrations.AlterField(
            model_name='siteobservation',
            name='ligand_mol_file',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='siteobservation',
            name='trans_matrix_info',
            field=models.JSONField(encoder=django.core.serializers.json.DjangoJSONEncoder, null=True),
        ),
        migrations.CreateModel(
            name='HistoricalSiteObservation',
            fields=[
                ('id', models.IntegerField(auto_created=True, blank=True, db_index=True, verbose_name='ID')),
                ('code', models.TextField()),
                ('bound_file', models.TextField(max_length=255, null=True)),
                ('apo_solv_file', models.TextField(max_length=255, null=True)),
                ('apo_desolv_file', models.TextField(max_length=255, null=True)),
                ('apo_file', models.TextField(max_length=255, null=True)),
                ('xmap_2fofc_file', models.TextField(max_length=255, null=True)),
                ('event_file', models.TextField(max_length=255, null=True)),
                ('artefacts_file', models.TextField(max_length=255, null=True)),
                ('trans_matrix_info', models.JSONField(encoder=django.core.serializers.json.DjangoJSONEncoder, null=True)),
                ('pdb_header_file', models.TextField(max_length=255, null=True)),
                ('smiles', models.TextField()),
                ('seq_id', models.IntegerField()),
                ('chain_id', models.CharField(max_length=1)),
                ('ligand_mol_file', models.TextField(null=True)),
                ('history_id', models.AutoField(primary_key=True, serialize=False)),
                ('history_date', models.DateTimeField()),
                ('history_change_reason', models.CharField(max_length=100, null=True)),
                ('history_type', models.CharField(choices=[('+', 'Created'), ('~', 'Changed'), ('-', 'Deleted')], max_length=1)),
                ('canon_site_conf', models.ForeignKey(blank=True, db_constraint=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='+', to='viewer.canonsiteconf')),
                ('cmpd', models.ForeignKey(blank=True, db_constraint=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='+', to='viewer.compound')),
                ('experiment', models.ForeignKey(blank=True, db_constraint=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='+', to='viewer.experiment')),
                ('history_user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='+', to=settings.AUTH_USER_MODEL)),
                ('xtalform_site', models.ForeignKey(blank=True, db_constraint=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='+', to='viewer.xtalformsite')),
            ],
            options={
                'verbose_name': 'historical site observation',
                'ordering': ('-history_date', '-history_id'),
                'get_latest_by': 'history_date',
            },
            bases=(simple_history.models.HistoricalChanges, models.Model),
        ),
    ]