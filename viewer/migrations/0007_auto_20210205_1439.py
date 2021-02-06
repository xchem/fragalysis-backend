# Generated by Django 3.1 on 2021-02-05 14:39

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0006_discoursecategory_discoursetopic'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='computedmolecule',
            name='pdb_info',
        ),
        migrations.AddField(
            model_name='computedmolecule',
            name='pdb',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.PROTECT, to='viewer.protein'),
        ),
    ]