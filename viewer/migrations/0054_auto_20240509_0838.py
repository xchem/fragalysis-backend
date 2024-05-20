# Generated by Django 3.2.24 on 2024-05-09 08:38

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0053_auto_20240503_1302'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='experiment',
            name='version',
        ),
        migrations.AddField(
            model_name='experimentupload',
            name='upload_data_dir',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='experimentupload',
            name='upload_version',
            field=models.PositiveSmallIntegerField(default=1),
        ),
    ]