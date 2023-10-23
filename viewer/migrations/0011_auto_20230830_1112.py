# Generated by Django 3.2.20 on 2023-08-30 11:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0010_remove_xtalform_xtalform_id'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='canonsite',
            name='quat_assembly',
        ),
        migrations.AddConstraint(
            model_name='xtalformquatassembly',
            constraint=models.UniqueConstraint(fields=('xtalform', 'quat_assembly', 'assembly_id'), name='unique_xtalformquatassembly'),
        ),
    ]