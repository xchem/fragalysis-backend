# Generated by Django 3.1 on 2023-02-01 19:15

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('car', '0004_auto_20221125_0856'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='reaction',
            name='recipetype',
        ),
        migrations.AddField(
            model_name='reaction',
            name='recipe',
            field=models.CharField(default='standard', max_length=50),
        ),
    ]
