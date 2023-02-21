# Generated by Django 3.1.14 on 2022-12-07 14:21

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('viewer', '0024_add_job_request_start_and_finish_times'),
    ]

    operations = [
        migrations.CreateModel(
            name='Squonk2Org',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('uuid', models.TextField(max_length=40)),
                ('name', models.TextField(max_length=80)),
                ('as_url', models.URLField()),
                ('as_version', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='Squonk2Unit',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('uuid', models.TextField(max_length=41)),
                ('name', models.TextField()),
                ('organisation', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.squonk2org')),
            ],
        ),
        migrations.CreateModel(
            name='Squonk2Project',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('uuid', models.TextField(max_length=44)),
                ('name', models.TextField()),
                ('product_uuid', models.TextField(max_length=44)),
                ('unit', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='viewer.squonk2unit')),
            ],
        ),
    ]