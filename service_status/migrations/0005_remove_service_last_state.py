# Generated by Django 3.2.25 on 2024-05-17 15:13

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('service_status', '0004_auto_20240517_1510'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='service',
            name='last_state',
        ),
    ]
