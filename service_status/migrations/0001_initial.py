# Generated by Django 3.2.25 on 2024-05-17 13:54

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='ServiceState',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('state', models.TextField()),
                ('display_name', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='Services',
            fields=[
                ('service', models.TextField(primary_key=True, serialize=False)),
                ('display_name', models.TextField()),
                ('frequency', models.PositiveSmallIntegerField(default=30, help_text='Ping frequency in seconds')),
                ('last_states_of_same_type', models.IntegerField(null=True)),
                ('last_query_time', models.DateTimeField(null=True)),
                ('last_success', models.DateTimeField(null=True)),
                ('last_failure', models.DateTimeField(null=True)),
                ('last_state', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='service_status.servicestate')),
            ],
        ),
    ]
