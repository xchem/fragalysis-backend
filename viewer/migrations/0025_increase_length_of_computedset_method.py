# Generated by Django 3.2.23 on 2023-12-19 14:31

import shortuuid.django_fields
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ('viewer', '0024_modified_computedset'),
    ]

    operations = [
        migrations.AlterField(
            model_name='computedset',
            name='method',
            field=models.CharField(
                blank=True,
                help_text='The name of the algorithmic method used to generate the compounds (e.g. Fragmenstein)',
                max_length=50,
                null=True,
            ),
        ),
        migrations.AlterField(
            model_name='computedmolecule',
            name='identifier',
            field=shortuuid.django_fields.ShortUUIDField(
                alphabet='ACDEFGHJKLMNPRSTUVWXYZ345679',
                blank=True,
                help_text="A four character string of non-confusing uppercase letters and digits for easy reference. This is combined with the Target to form the ComputedMolecule's name",
                length=4,
                max_length=4,
                null=True,
                prefix='',
            ),
        ),
    ]
