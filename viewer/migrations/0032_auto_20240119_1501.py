# Generated by Django 3.2.23 on 2024-01-19 15:01

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0031_compound_compound_code'),
    ]

    operations = [
        migrations.AddField(
            model_name='experimentupload',
            name='conformer_site_transforms',
            field=models.FileField(
                default='trans.yaml',
                max_length=255,
                upload_to='experiment-upload/',
            ),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='experimentupload',
            name='neighbourhood_transforms',
            field=models.FileField(
                default='trans.yaml', max_length=255, upload_to='experiment-upload/'
            ),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='experimentupload',
            name='reference_structure_transforms',
            field=models.FileField(
                default='trans.yaml', max_length=255, upload_to='experiment-upload/'
            ),
            preserve_default=False,
        ),
    ]
