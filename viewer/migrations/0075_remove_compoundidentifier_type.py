# Generated by Django 3.2.25 on 2024-10-18 10:38

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0074_auto_20241018_1028'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compoundidentifier',
            name='type',
        ),
    ]
