from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0156_userrole'),
    ]

    operations = [
        migrations.AddField(
            model_name='downloadlinks',
            name='expiry_reason',
            field=models.TextField(
                help_text='Why the download was expired (e.g. the task was lost).'
                ' Shown to users querying a failed download.',
                null=True,
            ),
        ),
    ]
