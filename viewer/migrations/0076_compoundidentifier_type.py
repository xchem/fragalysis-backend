# Generated by Django 3.2.25 on 2024-10-18 10:39

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0075_remove_compoundidentifier_type'),
    ]

    operations = [
        migrations.AddField(
            model_name='compoundidentifier',
            name='type',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='viewer.compoundidentifiertype'),
        ),
    ]