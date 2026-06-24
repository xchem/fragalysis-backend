# This is a mitigation for the issue (1609, 1921) where certain media
# files simply disappear. So far only .pdb, .cif and .mtz files from
# crystallographic_files/../ seem to be affected (all referenced in
# Experiment model. Reason for disappearance is not yet known, this
# command restores the files from the uploaded tarball.

import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

from django.conf import settings
from django.core.management.base import BaseCommand

from viewer.models import Target

EXP_FILES = ('pdb_info', 'mtz_info', 'cif_info')


class Command(BaseCommand):
    help = 'Restore files from v2 tarballs for targets converted to v3'

    dry_run = False

    def add_arguments(self, parser):
        parser.add_argument(
            '-d',
            '--dry-run',
            action="store_true",
            help="Don't restore, just list the files",
        )

    def handle(self, *args, **kwargs):
        # Unused args
        del args
        dry_run = kwargs.get('dry_run', False)

        for target in Target.objects.all():
            for exp_upload in target.experimentupload_set.all():
                tarball_path = Path(
                    settings.MEDIA_ROOT,
                    settings.TARGET_LOADER_MEDIA_DIRECTORY,
                    target.zip_archive.name,
                    exp_upload.file.name,
                )
                # Use the default temporary location (an emptyDir mounted at
                # /tmp within the Pod) rather than the shared MEDIA_ROOT volume;
                # restored files are copied into MEDIA_ROOT with shutil.copy().
                # See ticket #935.
                with TemporaryDirectory() as decompress_dir:
                    process = subprocess.Popen(
                        [
                            "tar",
                            "-I",
                            "pigz",
                            "-xf",
                            str(tarball_path),
                            "-C",
                            decompress_dir,
                        ],
                        stdout=subprocess.PIPE,
                    )
                    process.wait()
                    for exp in exp_upload.experiment_set.all():
                        # Experiment files are not altloc'ed, restore as they are

                        source_root = Path(
                            decompress_dir,
                            exp_upload.upload_data_dir,
                            'crystallographic_files',
                        )

                        for field in EXP_FILES:
                            model_attr = getattr(exp, field)
                            if model_attr and model_attr != 'None':
                                path = Path(settings.MEDIA_ROOT, model_attr.name)
                                if not path.exists():
                                    # find the file from the decomppresed dir
                                    source_path = source_root.joinpath(
                                        *Path(model_attr.name).parts[4:]
                                    )
                                    if source_path.exists():
                                        self.stdout.write(f'Restoring {path}')
                                        if not dry_run:
                                            shutil.copy(source_path, path)
                                    else:
                                        self.stdout.write(
                                            self.style.ERROR(
                                                f'Source path {source_path} does not exist'
                                            )
                                        )

                            # array field, special case
                            # although.. I don't recall these ever going missing
                            model_attr = getattr(exp, 'map_info')
                            if model_attr:
                                for f in model_attr:
                                    path = Path(settings.MEDIA_ROOT, f)
                                    if not path.exists():
                                        # find the file from the decomppresed dir
                                        source_path = source_root.joinpath(
                                            *Path(model_attr.name).parts[4:]
                                        )
                                        if source_path.exists():
                                            # these need to be copied
                                            # annoyingly, python
                                            # 3.14's pathlib can copy,
                                            # but we're not there yet
                                            self.stdout.write(f'Restoring {path}')
                                            if not dry_run:
                                                shutil.copy(source_path, path)

                                        else:
                                            self.stdout.write(
                                                self.style.ERROR(
                                                    f'Source path {source_path} does not exist'
                                                )
                                            )

                        # don't need site observation file restore,
                        # according to current stats, they're not
                        # going missing
