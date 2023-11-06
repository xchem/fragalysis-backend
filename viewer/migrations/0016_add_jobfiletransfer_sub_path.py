# Generated by Django 3.1.14 on 2023-11-06 14:23

from django.db import migrations
import shortuuid.django_fields


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0015_auto_20231023_1049'),
    ]

    operations = [
        migrations.AddField(
            model_name='jobfiletransfer',
            name='sub_path',
            field=shortuuid.django_fields.ShortUUIDField(alphabet='abcdefghijklmnopqrstuvwxyz', length=4, max_length=4, null=True, prefix=''),
        ),
    ]
