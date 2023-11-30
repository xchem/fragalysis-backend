# Generated by Django 3.2.20 on 2023-09-08 10:11

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0011_auto_20230830_1112'),
    ]

    operations = [
        migrations.AlterField(
            model_name='experimentupload',
            name='task_id',
            field=models.CharField(
                help_text='Celery task ID responsible for the upload',
                max_length=50,
                null=True,
            ),
        ),
    ]
