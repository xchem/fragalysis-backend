# Generated by Django 3.2.24 on 2024-05-03 13:02

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0052_auto_20240503_1255'),
    ]

    operations = [
        migrations.AddField(
            model_name='canonsite',
            name='version',
            field=models.PositiveSmallIntegerField(default=1),
        ),
        migrations.AddField(
            model_name='canonsiteconf',
            name='version',
            field=models.PositiveSmallIntegerField(default=1),
        ),
        migrations.AddField(
            model_name='xtalformsite',
            name='version',
            field=models.PositiveSmallIntegerField(default=1),
        ),
    ]
