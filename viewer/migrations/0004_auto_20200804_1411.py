# Generated by Django 3.1 on 2020-08-04 14:11

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0003_auto_20200804_1241'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='activitypoint',
            options={},
        ),
        migrations.AlterModelOptions(
            name='compound',
            options={},
        ),
        migrations.AlterModelOptions(
            name='molecule',
            options={},
        ),
        migrations.AlterModelOptions(
            name='project',
            options={},
        ),
        migrations.AlterModelOptions(
            name='protein',
            options={},
        ),
        migrations.AlterModelOptions(
            name='sessionproject',
            options={},
        ),
        migrations.AlterModelOptions(
            name='snapshot',
            options={'managed': True},
        ),
        migrations.AlterModelOptions(
            name='target',
            options={},
        ),
    ]