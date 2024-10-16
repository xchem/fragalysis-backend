# Generated by Django 3.2.25 on 2024-10-10 14:37

from pathlib import Path

from django.db import migrations

from fragalysis.settings import MEDIA_ROOT, TARGET_LOADER_MEDIA_DIRECTORY

from viewer.utils import sanitize_directory_name


# continuation of the data migration to replace m2m between target and
# project with fk in target to project
# starts in 0066, where target gets the fk to project populated
# migrations in between take care of changing the schema
# this one now moves the data to correct locations

class Migration(migrations.Migration):

    dependencies = [
        ('viewer', '0071_target_unique_target_in_project'),
    ]


    def move_data(apps, schema_editor):
        Target = apps.get_model('viewer', 'Target')
        root_path = Path(MEDIA_ROOT).joinpath(TARGET_LOADER_MEDIA_DIRECTORY)
        for target in Target.objects.all():


            if target.zip_archive:
                old_data_path = root_path.joinpath(target.zip_archive.name)
                new_archive_dir = sanitize_directory_name(
                    f"{target.zip_archive.name}_{target.project.title}",
                    root_path,
                )
                new_data_path = root_path.joinpath(new_archive_dir)
                old_data_path.rename(new_data_path)
                # move tarballs
                for exp_upload in target.experimentupload_set.all():
                    tarball_path = root_path.joinpath(exp_upload.file.name)
                    new_tarball_path = new_data_path.joinpath(exp_upload.file.name)
                    tarball_path.rename(new_tarball_path)
            else:
                # old system, data path contains task_id
                new_archive_dir = sanitize_directory_name(
                    f"{target.title}_{target.project.title}",
                    root_path,
                )
                new_data_path = root_path.joinpath(new_archive_dir)
                new_data_path.mkdir(parents=True, exist_ok=False)
                for exp_upload in target.experimentupload_set.all():
                    old_data_path = root_path.joinpath(str(exp_upload.task_id))
                    # can only be one
                    target_subdir = [p for p in old_data_path.glob("*") if p.is_dir()][0]
                    target_tarball = [p for p in old_data_path.glob("*") if p.is_file()][0]

                    upload_dir = next(target_subdir.glob('upload_*'))
                    # NB! move tarball first!
                    new_tarball_path = new_data_path.joinpath(target_tarball.name)
                    target_tarball.rename(new_tarball_path)
                    # and then the rest of the data
                    upload_dir.rename(new_data_path.joinpath(upload_dir.parts[-1]))



            target.zip_archive = new_archive_dir
            target.save()


    def unmove_data(apps, schema_editor):
        pass

    operations = [
        migrations.RunPython(move_data, unmove_data)
    ]
