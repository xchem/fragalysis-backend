"""
download_structures.py

Methods for downloading a Target Zip file used by the download_structure API.
"""

import copy
import csv
import json
import logging
import os
import shutil
import subprocess
import time
import zipfile
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from enum import Enum
from io import BytesIO, StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict

import humanize
import pandas as pd
import pandoc
import requests
from celery.result import AsyncResult
from django.conf import settings
from django.db.models import Exists, F, OuterRef, Value
from django.db.models.fields import CharField
from django.db.models.functions import Concat
from rdkit import Chem

from api.infections import INFECTION_STRUCTURE_DOWNLOAD_TASK, have_infection
from fragalysis.celery import app as celery_app
from viewer.models import DownloadLinks, SiteObservation
from viewer.utils import clean_filename

from .logger_adapters import TaskLoggerAdapter
from .tags import get_metadata_fields

logger = logging.getLogger(__name__)

# Length of time to keep records of dynamic links before "expiry"
KEEP_UNTIL_DURATION = timedelta(minutes=settings.DOWNLOAD_KEEP_UNTIL_DURATION_M)
# Length of time to "expired" records
HARD_EXPIRY_GRACE_PERIOD = timedelta(minutes=settings.HARD_EXPIRY_GRACE_PERIOD_M)
# An in-progress download is given this long to start before the "lost task"
# housekeeping pass will consider it (see expire_lost_download_records).
TASK_START_GRACE_PERIOD = timedelta(minutes=settings.DOWNLOAD_TASK_START_GRACE_M)
# The longest an in-progress download may run before it is presumed lost (e.g.
# the worker was restarted, leaving the Celery result frozen "in progress").
TASK_MAX_RUNTIME = timedelta(minutes=settings.DOWNLOAD_TASK_MAX_RUNTIME_M)

# Filepaths mapping for writing associated files to the zip archive.
# Note that if this is set to 'aligned' then the files will be placed in
# the protein code subdirectory of the aligned directory
# (as for the target upload).
_ZIP_FILEPATHS = {
    'apo_file': ('aligned'),  # SiteObservation: apo_file
    'apo_solv_file': ('aligned'),  # SiteObservation: apo_solv_file
    'apo_desolv_file': ('aligned'),  # SiteObservation: apo_desolv_file
    'bound_file': ('aligned'),  # SiteObservation: bound_file
    'sdf_info': ('aligned'),  # SiteObservation: ligand_mol_file (indirectly)
    'ligand_sdf': ('aligned'),  # SiteObservation: ligand_sdf
    'ligand_smiles': ('aligned'),  # SiteObservation: ligand_smiles
    'ligand_pdb': ('aligned'),  # SiteObservation: ligand_pdb
    'smiles_info': (''),  # SiteObservation: smiles_info (indirectly)
    # those above are all controlled by serializer's all_aligned_structures flag
    'sigmaa_file': ('aligned'),  # SiteObservation: sigmaa_file
    'diff_file': ('aligned'),  # SiteObservation: diff_file
    'event_file': ('aligned'),  # SiteObservation: ligand_pdb
    'pdb_info': ('aligned'),  # Experiment: cif_info
    'cif_info': ('aligned'),  # Experiment: cif_info
    'mtz_info': ('aligned'),  # Experiment: mtz_info
    'map_info': ('aligned'),  # Experiment: map_info (multiple files)
    'single_sdf_file': (''),
    'metadata_info': (''),
    'trans_matrix_info': (''),
    'extra_files': ('extra_files'),
    'readme': (''),
}

# {urls:path in archive} to scripts to be automatically included in downloads
_SCRIPTS = {
    'https://github.com/xchem/fragalysis-pymol-scripts': 'pymol',
}


@dataclass(frozen=True, eq=True)
class ArchiveFile:
    path: str
    archive_path: str
    site_observation: SiteObservation | None = None


class ProcessState(str, Enum):
    """Downloader progress state.

    PROCESSING - all good, upload in progress
    SUCCESS    - processing complete, all good
    FAILED     - processing complete, failed
    """

    PROCESSING = "PROCESSING"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"


# Dictionary containing all references needed to create the zip file
# NB you may need to add a version number to this at some point...
zip_template = {
    'proteins': {
        'apo_file': {},
        'apo_solv_file': {},
        'apo_desolv_file': {},
        'bound_file': {},
        'pdb_info': {},
        'cif_info': {},
        'mtz_info': {},
        'map_info': {},
        'event_file': {},
        'diff_file': {},
        'sigmaa_file': {},
        'ligand_pdb': {},
        'ligand_sdf': {},
        'ligand_smiles': {},
        # additional ccp4 files, issue 1448
        'event_file_crystallographic': {},
        'diff_file_crystallographic': {},
        'sigmaa_file_crystallographic': {},
    },
    'molecules': {
        'sdf_files': {},
        'sdf_info': False,
        'single_sdf_file': False,
        'smiles_info': {},
    },
    'metadata_info': None,
    'trans_matrix_info': None,
    'compound_sets': None,
    'soakdb_files': None,
}


_ERROR_FILE = 'errors.csv'

# unlike v1, metadata doesn't exist anymore, needs compiling
_METADATA_FILE = 'metadata.csv'


# NB! I believe the download archive creation can be massively
# simplified by using database operations. Most of the fields can be
# added to SiteObservation queryset with annotations. There is a
# map_info field which currently is ArrayField which may be a bit of a
# problem. This could be normalized or perhaps using Pefetch objects
# below. And then, potentially, the same queryset re-used for
# metadata. This is a bit of work and was out of scope for this ticket
# (1854) but potentially worth revisiting.

# site_observations = site_observations.annotate(
#     # this is the key to use later for name substitutions
#     longlongcode=Concat(
#         F('experiment__code'),
#         Value('_'),
#         F('chain_id'),
#         Value('_'),
#         F('seq_id'),
#         Value('_'),
#         F('version'),
#         Value('_'),
#         F('canon_site_conf__canon_site__name'),
#         Value('+'),
#         F('canon_site_conf__canon_site__version'),
#         output_field=CharField(),
#     ),
# )

# prefetch = Prefetch(
#     "cmpd__all_identifiers",
#     queryset=CompoundIdentifier.objects.annotate(
#         substituted_name=Replace(
#             F("name"),
#             Value("the_"),
#             Value("a___"),
#             output_field=CharField(),
#         ),
#     ),
#     to_attr="annotated_name"
# )


class DownloadStructures:
    def __init__(
        self,
        *,
        tempdir,
        task,
        target,
        use_zip,
        target_access_string,
    ):
        self.task = task
        self.use_zip = use_zip
        self.tempdir = tempdir
        self.target_access_string = target_access_string
        self.target = target

        self._temp_path = Path(self.tempdir)
        self._combined_sdf_path = self._temp_path.joinpath(
            f'{target.title}_combined.sdf'
        )
        self._error_file = self._temp_path.joinpath(_ERROR_FILE)
        # A custom logging adapter to refine the standard logger
        # by adding lightweight context to the messages we log.
        # We don't add Target, TAS, or Username.
        #
        # The assumption is that it's already logged by the outer
        # 'create_download_link()' function and probably unnecessary to add this
        # information on every line - the Task ID should be enough to find matching
        # lines.
        self._logger = TaskLoggerAdapter(
            logger,
            {
                'task': str(task.request.id),
                'marker': 'DOWNLOAD',
            },
        )

    @property
    def temp_path(self) -> Path:
        return self._temp_path

    @property
    def combined_sdf_path(self) -> Path:
        return self._combined_sdf_path

    @property
    def error_file(self) -> Path:
        return self._error_file

    def create_content_dict(
        self,
        *,
        site_obvs,
        protein_params,
        other_params,
    ):
        """Create dict containing the listing of the tarball contents.

        Args:
            target
            proteins
            protein_params
            other_params

        Returns:
            [dict]: [dictionary containing the file contents]
        """
        self._logger.info('Processing %d SiteObservations', site_obvs.count())
        self.update_task(ProcessState.PROCESSING, 'Creating tarball contents...')

        # Read through zip_params to compile the parameters
        zip_contents: Dict[str, Any] = copy.deepcopy(zip_template)
        site_obvs = site_obvs.annotate(
            # would there be any point in
            # a) adding a method to SiteObservation model_attr
            # b) adding the value to database directly?
            longlongcode=Concat(
                F('experiment__code'),
                Value('_'),
                F('chain_id'),
                Value('_'),
                F('seq_id'),
                Value('_'),
                F('altloc'),
                Value('_'),
                F('version'),
                Value('_'),
                F('canon_site_conf__canon_site__name'),
                Value('+'),
                F('canon_site_conf__canon_site__version'),
                output_field=CharField(),
            ),
        )
        for so in site_obvs:
            for param in protein_params:
                if protein_params[param] is True:
                    if param in ['pdb_info', 'mtz_info', 'cif_info', 'map_info']:
                        # experiment object
                        model_attr = getattr(so.experiment, param)
                        self._logger.debug(
                            'Adding param to zip: %s, value: %s', param, model_attr
                        )
                        if param != 'map_info':
                            # treat all params as list
                            model_attr = (
                                [model_attr.name]
                                # None - some weird glitch in storing the values
                                if model_attr and not str(model_attr).find('None') > -1
                                else [param]
                            )

                        afile = []
                        for f in model_attr:
                            apath = Path('crystallographic_files', so.experiment.code)
                            if model_attr and model_attr != 'None':
                                archive_path = str(apath.joinpath(Path(f).parts[-1]))
                            else:
                                archive_path = str(apath.joinpath(param))
                            afile.append(ArchiveFile(path=f, archive_path=archive_path))

                        # changed in [can't find the ticket
                        # number]. previously zip_contents was updated
                        # in the end of for cycle instead of
                        # crystallographic and aligned separately, but
                        # in attempt to deduplicate files in
                        # crystallographic dir, I'm moving the update
                        # statement to each category's section. I'm a
                        # bit worried it might silently overwrite some
                        # files but I can't think of any other than
                        # those I actually want to
                        zip_contents['proteins'][param][so.experiment.code] = afile

                    elif param in [
                        'bound_file',
                        'apo_file',
                        'apo_solv_file',
                        'apo_desolv_file',
                        'sigmaa_file',
                        'event_file',
                        'artefacts_file',
                        'pdb_header_file',
                        'ligand_pdb',
                        'ligand_sdf',
                        'ligand_smiles',
                        'diff_file',
                    ]:
                        # siteobservation object

                        model_attr = getattr(so, param)
                        self._logger.debug(
                            'Adding param to zip: %s, value: %s', param, model_attr
                        )
                        apath = Path('aligned_files').joinpath(so.code)
                        if model_attr and model_attr != 'None':
                            archive_path = str(
                                apath.joinpath(
                                    Path(model_attr.name)
                                    .parts[-1]
                                    .replace(so.longlongcode, so.code)
                                )
                            )
                        else:
                            # file not in upload
                            archive_path = str(apath.joinpath(param))

                        afile = [
                            ArchiveFile(
                                path=model_attr.name,
                                archive_path=archive_path,
                            )
                        ]

                        zip_contents['proteins'][param][so.code] = afile

                    else:
                        self._logger.warning('Unexpected param: %s', param)
                        continue

                    # add additional ccp4 files (issue 1448)
                    ccps = ('sigmaa_file', 'diff_file', 'event_file')
                    if param in ccps:
                        # these only come from siteobservation object
                        model_attr = getattr(so, param)
                        if model_attr and model_attr != 'None':
                            apath = Path('aligned_files').joinpath(so.code)
                            ccp_path = Path(model_attr.name)
                            path = ccp_path.parent.joinpath(
                                f'{ccp_path.stem}_crystallographic{ccp_path.suffix}'
                            )
                            archive_path = str(
                                apath.joinpath(
                                    path.parts[-1].replace(so.longlongcode, so.code)
                                )
                            )

                            afile = [
                                ArchiveFile(
                                    path=str(path),
                                    archive_path=archive_path,
                                )
                            ]
                            zip_contents['proteins'][f'{param}_crystallographic'][
                                so.code
                            ] = afile

        zip_contents['molecules']['single_sdf_file'] = other_params['single_sdf_file']
        zip_contents['molecules']['sdf_info'] = other_params['sdf_info']

        if other_params['sdf_info'] or other_params['single_sdf_file']:
            num_molecules_collected = 0
            num_missing_sd_files = 0
            for so in site_obvs:
                if so.ligand_sdf:
                    # There is an SD file (normal)
                    archive_path = str(
                        Path('aligned_files')
                        .joinpath(so.code)
                        .joinpath(f'{so.code}.sdf')
                    )
                    file_path = str(
                        Path(settings.MEDIA_ROOT).joinpath(so.ligand_sdf.name)
                    )
                    # mypy, you make me do stupid things..
                    if Path(file_path).exists():
                        # expects str
                        file_path = str(file_path)
                    else:
                        file_path = ''
                    # path is ignored when writing sdfs but mandatory field
                    zip_contents['molecules']['sdf_files'].update(
                        {
                            ArchiveFile(
                                path=file_path,
                                archive_path=archive_path,
                                site_observation=so,
                            ): so.code
                        }
                    )
                    num_molecules_collected += 1
                else:
                    # No file value (odd).
                    self._logger.warning(
                        "SiteObservation record's 'ligand_sdf' isn't set (%s)", so
                    )
                    num_missing_sd_files += 1

            # Report (in the log) anomalies
            if num_molecules_collected == 0:
                self._logger.warning('No SD files collected')
            else:
                self._logger.info('%s SD files collected', num_molecules_collected)

            if site_obvs.count() != num_molecules_collected:
                self._logger.warning(
                    'Expected %d files, got %d',
                    site_obvs.count(),
                    num_molecules_collected,
                )

            if num_missing_sd_files > 0:
                self._logger.warning('%d missing files', num_missing_sd_files)

        # The smiles at molecule level may not be unique.
        if other_params['smiles_info'] is True:
            for molecule in site_obvs:
                zip_contents['molecules']['smiles_info'].update({molecule.smiles: None})

        zip_contents['metadata_info'] = other_params['metadata_info']

        # Add the trans matrix files
        zip_contents['trans_matrix_info'] = other_params['trans_matrix_info']
        zip_contents['compound_sets'] = other_params['compound_sets']
        zip_contents['soakdb_files'] = other_params['soakdb_files']

        return zip_contents

    def create_tarball(
        self,
        *,
        zip_contents,
        file_url,
        original_search,
        site_observations,
    ):
        """Write a ZIP file containing data from an input dictionary."""

        self._logger.info('file_url="%s"', file_url)
        self._logger.info(
            'single_sdf_file="%s"', zip_contents['molecules']['single_sdf_file']
        )
        self._logger.debug('sdf_files=%s', zip_contents['molecules']['sdf_files'])

        self._logger.debug('zip_contents=%s', zip_contents)
        self.update_task(ProcessState.PROCESSING, 'Creating tarball...')

        download_path = os.path.dirname(file_url)
        self._logger.info('Creating download path (%s)', download_path)
        os.makedirs(download_path, exist_ok=True)

        error_filename = str(self.error_file)
        error_file = open(error_filename, "w", encoding="utf-8")
        error_file.write("Param,Code,File not found when assembling download\n")
        errors = 0

        # If a single sdf file is also wanted then create file to
        # add sdf files to a file called {target}_combined.sdf.
        combined_sdf_file = None
        if zip_contents['molecules']['single_sdf_file'] is True:
            combined_sdf_file = str(self.combined_sdf_path)
            self._logger.info('combined_sdf_file=%s', combined_sdf_file)

        # Read through zip_contents to compile the file
        self.update_task(ProcessState.PROCESSING, 'Adding PDBs...')
        errors += self._protein_files_zip(zip_contents, error_file)
        if errors > 0:
            self._logger.warning('After _protein_files_zip() errors=%s', errors)

        self.update_task(ProcessState.PROCESSING, 'Adding SDFs...')
        if zip_contents['molecules']['sdf_files']:
            errors_before = errors
            errors += self._molecule_files_zip(
                zip_contents, combined_sdf_file, error_file
            )
            if errors > errors_before:
                self._logger.warning('After _molecule_files_zip() errors=%s', errors)

        # If smiles info is required, then write one column for each molecule
        # to a smiles.smi file and then add to the archive.
        if zip_contents['molecules']['smiles_info']:
            self.update_task(ProcessState.PROCESSING, 'Adding SMILES...')
            self._smiles_files_zip(zip_contents)

        # compile and add metadata.csv
        if zip_contents['metadata_info']:
            self.update_task(ProcessState.PROCESSING, 'Adding metadata.csv...')
            self._metadata_file_zip(self.target, site_observations)

        if zip_contents['trans_matrix_info']:
            self.update_task(
                ProcessState.PROCESSING, 'Adding transformation matrix files...'
            )
            self._trans_matrix_files_zip(self.target)

        self.update_task(ProcessState.PROCESSING, 'Adding extra files...')
        self._extra_files_zip(self.target, soakdb_files=zip_contents['soakdb_files'])

        self.update_task(ProcessState.PROCESSING, 'Adding YAMLs...')
        self._yaml_files_zip(
            self.target, transforms_requested=zip_contents['trans_matrix_info']
        )

        self.update_task(ProcessState.PROCESSING, 'Adding scripts...')
        self._additional_scripts_zip(_SCRIPTS)

        self.update_task(ProcessState.PROCESSING, 'Adding compound sets...')
        if zip_contents['compound_sets']:
            self._compound_sets_zip(self.target)

        error_file.close()

        # memo to self: this function includes the file list in
        # the result, so needs to come last
        self.update_task(ProcessState.PROCESSING, 'Creating documentation...')
        self._document_file_zip(original_search)

        import timeit

        t_0 = timeit.default_timer()
        self.compress_directory(self.temp_path, file_url)
        t_end = timeit.default_timer()
        self._logger.debug('timings, compression time: %s', t_end - t_0)

    def update_task(self, status: ProcessState, message: str):
        self.task.update_state(
            state=status,
            meta={
                "proposal_ref": self.target_access_string,
                "description": message,
            },
        )

    def write_symlink(self, path, archive_name):
        archive_path = self.temp_path.joinpath(archive_name)
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        archive_path.symlink_to(Path(path))

    def write_file(self, contents, archive_name):
        archive_path = self.temp_path.joinpath(archive_name)
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        archive_path.write_text(str(contents), encoding='utf-8')

    def _additional_scripts_zip(self, scripts) -> None:
        for script_url, script_path in scripts.items():
            zip_url = script_url.rstrip('/') + '/archive/refs/heads/main.zip'
            response = requests.get(zip_url)
            try:
                response.raise_for_status()
            except Exception as exc:
                self.write_file(str(exc), f'scripts/ERROR_DOWNLOADING_{script_path}')

            try:
                with zipfile.ZipFile(BytesIO(response.content)) as zip_file:
                    script_path = Path('scripts', script_path)
                    for zip_info in zip_file.infolist():
                        if zip_info.is_dir():
                            continue  # Skip directories
                        file_data = zip_file.read(zip_info.filename)
                        archive_path = script_path.joinpath(
                            *Path(zip_info.filename).parts[1:]
                        )
                        self.write_file(file_data, str(archive_path))
            except zipfile.BadZipFile as exc:
                self.write_file(str(exc), f'scripts/ERROR_DOWNLOADING_{script_path}')

    def _add_empty_file(self, archive_path):
        """When file is missing, add an empty file to the archive.


        Used to send an explicit signal to the downloader that the file is
        missing.
        """
        self._logger.debug('+_add_empty_file: %s', archive_path)
        self.write_file('', f'{archive_path}_FILE_NOT_IN_UPLOAD')

    def _add_file_to_zip_aligned(self, code, archive_file):
        """Add the requested file to the zip archive.

        If the file is an SDF or MOL we insert the name of the molecule
        if it is not set  (the basename of the file).

        Args:
            ziparchive: Handle of zip archive
            code: original protein code stripped of any alternate name.
            filepath: filepath from record

        Returns:
            [boolean]: [True of record added to archive]
        """
        if not archive_file:
            # Odd - assume success
            self._logger.debug('No archive_file')
            return True

        # calling str on archive_file.path because could be None
        filepath = str(Path(settings.MEDIA_ROOT).joinpath(str(archive_file.path)))
        self._logger.debug(
            'value and type of archive path: %s, %s',
            archive_file.path,
            type(archive_file.path),
        )
        if archive_file.path:
            if Path(filepath).is_file():
                if _is_mol_or_sdf(filepath):
                    # It's a MOL or SD file.
                    # Read and (potentially) adjust the file
                    # and add to the archive as a string.
                    content = _read_and_patch_molecule_name(
                        filepath, molecule_name=code
                    )
                    self.write_file(content, archive_file.archive_path)
                else:
                    # Copy the file without modification
                    self.write_symlink(filepath, archive_file.archive_path)
                return True
            elif archive_file.site_observation:
                self.write_file(
                    _read_and_patch_molecule_name(
                        filepath, archive_file.site_observation
                    ),
                    archive_file.archive_path,
                )
                return True

        self._logger.debug('filepath "%s" is not a file', filepath)
        self._add_empty_file(archive_file.archive_path)

        return False

    def _add_file_to_sdf(self, combined_sdf_file, archive_file):
        """Append the requested sdf file to the single sdf file provided.

        Args:
            combined_sdf: Handle of combined_sdf_file
            filepath: filepath from record

        Returns:
            [boolean]: [True of record added]
        """
        if not archive_file.path:
            # Odd - assume success
            self._logger.debug('No archive_file.path value')
            return True

        if archive_file.path and archive_file.path != 'None':
            with open(combined_sdf_file, 'a', encoding='utf-8') as f_out:
                patched_sdf_content = _read_and_patch_molecule_name(
                    archive_file.path, archive_file.site_observation
                )
                f_out.write(patched_sdf_content)
            return True
        else:
            self._logger.debug('filepath "%s" is not a file', archive_file.path)

        return False

    def _protein_files_zip(self, zip_contents, error_file):
        """Write all protein related data to the ZIP file
        Returns protein errors
        """
        prot_errors = 0
        for param, files in zip_contents['proteins'].items():
            if not files:
                continue

            for prot, prot_file in files.items():
                for f in prot_file:
                    # memo to self: f is ArchiveFile object
                    if not self._add_file_to_zip_aligned(prot, f):
                        error_file.write(f'{param},{prot},{f.archive_path}\n')
                        prot_errors += 1

        return prot_errors

    def _molecule_files_zip(self, zip_contents, combined_sdf_file, error_file):
        """Write molecule (SD file) related data to the ZIP file
        Returns molecule errors
        """

        mol_errors = 0
        self._logger.info(
            'len(molecules.sd_files)=%s', len(zip_contents['molecules']['sdf_files'])
        )
        for archive_file, prot in zip_contents['molecules']['sdf_files'].items():
            # Do not try and process any missing SD files.
            if not archive_file:
                error_file.write(f'sdf_files,{prot},missing\n')
                mol_errors += 1
                continue

            if zip_contents['molecules'][
                'sdf_info'
            ] is True and not self._add_file_to_zip_aligned(
                prot.split(":")[0], archive_file
            ):
                error_file.write(f'sdf_info,{prot},{archive_file.path}\n')
                mol_errors += 1

            # Append sdf file on the Molecule record to the combined_sdf_file.
            if zip_contents['molecules'][
                'single_sdf_file'
            ] is True and not self._add_file_to_sdf(combined_sdf_file, archive_file):
                error_file.write(f'single_sdf_file,{prot},{archive_file.path}\n')
                mol_errors += 1

        return mol_errors

    def _smiles_files_zip(self, zip_contents):
        """Create and write the smiles file to the ZIP file"""
        smiles_filename = self.temp_path.joinpath('smiles.smi')
        self._logger.info('Creating SMILES file "%s"...', smiles_filename)

        num_smiles = 0
        with open(smiles_filename, 'w', encoding='utf-8') as smilesfile:
            for smi in zip_contents['molecules']['smiles_info']:
                self._logger.debug('Adding "%s"...', smi)
                smilesfile.write(f'{smi},')
                num_smiles += 1

        self._logger.info('Added %s SMILES', num_smiles)

    def _trans_matrix_files_zip(self, target):
        """Add transformation matrices to archive.

        Note that this will always be the latest information - even for
        preserved searches.
        """
        self._logger.info('+ Processing trans matrix files')

        # grab the last set of files for this target
        experiment_upload = target.experimentupload_set.order_by(
            'commit_datetime'
        ).last()

        trans_matrix_files = [
            f
            for f in (
                experiment_upload.neighbourhood_transforms,
                experiment_upload.conformer_site_transforms,
                experiment_upload.assembly_transforms,
                experiment_upload.reference_structure_transforms,
            )
            if f.name is not None
        ]
        for tmf in trans_matrix_files:
            filepath = Path(settings.MEDIA_ROOT).joinpath(str(tmf))
            archive_path = os.path.join(
                _ZIP_FILEPATHS['trans_matrix_info'],
                Path(str(tmf)).name,
            )
            if filepath.is_file():
                self.write_symlink(filepath, archive_path)
            else:
                self._logger.debug('File %s does not exist', Path(str(tmf)).name)
                self._add_empty_file(archive_path)

    def _metadata_file_zip(self, target, site_observations):
        """Compile and add metadata file to archive."""
        self._logger.info('Processing metadata...')

        header, annotations, values = get_metadata_fields(target)

        # fmt: off
        qs = SiteObservation.filter_manager.by_target(
            target=target,
        ).prefetch_related(
            'cmpd',
            # this wasn't a problem, until I tried to debug it and split
            # the qs (values later). why does this break it?
            # 'siteobservationtags',
        ).annotate(
            downloaded=Exists(
                site_observations.filter(
                    pk=OuterRef('pk'),
                ),
            )
        ).annotate(
            **annotations
        ).values(
            *values
        )
        # fmt: on

        df = pd.DataFrame(qs)
        self._logger.debug('qs: %s', qs)
        self._logger.debug('annotations: %s', annotations.keys())
        self._logger.debug('values: %s', values)

        columns = [header[values.index(k)] for k in df.columns]
        df.columns = columns

        buff = StringIO()
        df.to_csv(
            buff,
            header=True,
            index=False,
            encoding='utf-8',
            quoting=csv.QUOTE_NONNUMERIC,
            lineterminator="\n",
        )
        buff.seek(0)
        self.write_file(buff.getvalue(), _METADATA_FILE)

    def _extra_files_zip(self, target, soakdb_files=True):
        """If an extra info folder exists at the target root level, then
        copy the contents to the output file as is.
        Note that this will always be the latest information - even for
        preserved searches.
        """

        num_processed = 0
        num_extra_dir = 0
        # taking the latest upload for now

        experiment_upload = target.experimentupload_set.order_by(
            'commit_datetime'
        ).last()
        extra_files = (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(target.zip_archive.name)
            .joinpath(experiment_upload.upload_data_dir)
        )

        extra_files = extra_files.joinpath('extra_files')

        self._logger.debug('extra_files path 2: %s', extra_files)
        self._logger.info('Processing extra files (%s)...', extra_files)

        if extra_files.is_dir():
            num_extra_dir = num_extra_dir + 1

            for src_path in extra_files.rglob("*"):
                if src_path.is_file():
                    filepath = src_path.relative_to(extra_files)

                    if soakdb_files or (
                        not soakdb_files and str(filepath).find('soakdb_') < 0
                    ):
                        self._logger.info('Adding extra file "%s"...', src_path)
                        self.write_symlink(
                            src_path,
                            os.path.join(
                                f'{_ZIP_FILEPATHS["extra_files"]}_{num_extra_dir}',
                                filepath,
                            ),
                        )
                        num_processed += 1

        else:
            self._logger.info('Directory does not exist (%s)...', extra_files)

        if num_processed == 0:
            self._logger.info('No extra files found')
        else:
            self._logger.info('Processed %s extra files', num_processed)

    def _yaml_files_zip(self, target, transforms_requested: bool = False) -> None:
        """Add all yaml files (except transforms) from upload to ziparchive"""

        for experiment_upload in target.experimentupload_set.all():
            yaml_paths = (
                Path(settings.MEDIA_ROOT)
                .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
                .joinpath(target.zip_archive.name)
                .joinpath(experiment_upload.upload_data_dir)
            )

            transforms = [
                Path(f.name).name
                for f in (
                    experiment_upload.assembly_transforms,
                    experiment_upload.neighbourhood_transforms,
                    experiment_upload.conformer_site_transforms,
                    experiment_upload.reference_structure_transforms,
                )
                if f.name is not None
            ]

            archive_path = Path('yaml_files').joinpath(yaml_paths.parts[-1])

            yaml_files = [
                f
                for f in list(yaml_paths.glob("*.yaml"))
                if f.is_file() and f.name not in transforms
            ]

            self._logger.debug('Processing yaml files (%s)...', yaml_files)

            for file in yaml_files:
                self._logger.debug('Adding yaml file "%s"...', file)
                if not transforms_requested and file.name == 'neighbourhoods.yaml':
                    # don't add this file if transforms are not requested
                    continue
                self.write_symlink(file, str(Path(archive_path).joinpath(file.name)))

    def _compound_sets_zip(self, target) -> None:
        """Add compound sets to download"""

        self._logger.info('Processing computed sets...')
        for cset in target.computedset_set.all():
            archive_path = Path('virtual_hits').joinpath(cset.submitted_sdf.name)
            buff = StringIO()
            writer = Chem.SDWriter(buff)
            for cmol in cset.computed_molecules.all():
                self._logger.debug('Processing computed molecule (%s)...', cmol.name)
                mol = Chem.MolFromMolBlock(cmol.sdf_info)
                self._logger.debug('mol: %s', mol)
                mol.SetProp('_Name', cmol.name)
                for prop in cmol.numericalscorevalues_set.all():
                    mol.SetProp(prop.score.name, str(prop.value))
                for prop in cmol.textscorevalues_set.all():
                    mol.SetProp(prop.score.name, prop.value)

                writer.write(mol)

            self.write_file(buff.getvalue(), str(Path(archive_path)))

    def _document_file_zip(self, original_search):
        """Create the document file
        This consists of a template plus an added contents description.
        """

        self._logger.info('Creating documentation...')

        template_file = os.path.join(
            "/code/doc_templates", "download_readme_template.md"
        )
        readme_filepath = self.temp_path.joinpath('README.md')
        with open(str(readme_filepath), "a", encoding="utf-8") as readme:
            self._build_readme(readme, original_search, template_file)

        # Convert markdown to pdf file. Pandoc invokes external tools (the
        # pandoc binary, latex) that can fail for many reasons we don't want
        # to abort the whole download for — log and emit a placeholder PDF
        # carrying the exception text so the archive still contains README.pdf.
        pdf_filepath = self.temp_path.joinpath('README.pdf')
        try:
            doc = pandoc.read(open(readme_filepath, "r", encoding="utf-8").read())
            pandoc.write(
                doc, file=pdf_filepath, format='latex', options=["--columns=72"]
            )
        except Exception as exc:  # pylint: disable=broad-except
            self._logger.warning(
                'Pandoc README PDF generation failed (%s); '
                'writing placeholder README.pdf',
                exc,
            )
            _write_text_pdf(
                pdf_filepath,
                'README.pdf could not be generated.\n'
                '\n'
                'Pandoc raised the following error:\n'
                f'{exc}\n'
                '\n'
                'Please refer to README.md (in the same archive) for the '
                'full document.',
            )

        # self.write_symlink(pdf_filepath, os.path.join(_ZIP_FILEPATHS['readme'], 'README.pdf'))

    def _build_readme(self, readme, original_search, template_file):
        readme.write("# Documentation for the downloaded zipfile\n")
        # Download links
        readme.write("## Download details\n")
        # Removed as the URL wasn't being generated correctly.
        # readme.write("### Download URLs\n")
        # readme.write("- Download URL: <")
        # ext_url = _get_external_download_url(download_path, host)
        # readme.write(ext_url+">\n")

        # Original Search
        readme.write("\n### Download command (JSON)\n")
        readme.write(
            "JSON command sent from front-end to backend "
            "to generate the download. This can be reused "
            "programmatically as a POST command:\n\n"
        )
        readme.write(f"```{json.dumps(original_search)}" + "```\n\n")

        # Download Structure from the template
        # (but prepare for the template file not existing)?
        if os.path.isfile(template_file):
            with open(template_file, "r", encoding="utf-8") as template:
                readme.write(template.read())
        else:
            self._logger.warning('Could not find template file (%s)', template_file)

        # Files Included
        list_of_files = list(self.temp_path.rglob('*'))
        readme.write("\n## Files included\n")
        list_of_files.sort()
        for filename in list_of_files:
            readme.write(f'- {filename}' + '\n')

    def compress_directory(self, data_path, tarball_path):
        """Compress data for download.

        Two methods available, gzip (with pigz) and zip (7z).

        Compression levels were determined by running incomprehensive
        and non-conclusive local tests.
        """

        def check_popen_progress(compress_process, frequency):
            while compress_process.poll() is None:
                try:
                    current_size = os.path.getsize(tarball_path)
                except FileNotFoundError:
                    current_size = 0

                progress = min(current_size / directory_total_size_bytes, 1.0)
                self.update_task(
                    ProcessState.PROCESSING, f'Compressing tarball: {progress:.1%}'
                )
                time.sleep(frequency)

        directory_total_size_bytes = get_total_size(str(data_path))
        self._logger.info(
            'Creating tarball from %s (directory_total_size_bytes=%s)...',
            data_path,
            humanize.naturalsize(directory_total_size_bytes, binary=True),
        )
        poll_frequency = 2

        if self.use_zip:
            # NB! this is not the same tempdir where the symlinks are,
            # this is where the symlinks are resolved for 7z
            with TemporaryDirectory() as tmpdir:
                # unlike gzip, 7z cannot resolve symlinks. Have to do
                # this myself. Chose to use rsync instead of writing
                # files directly because
                # - would have to choose every time whether to use
                #   symlink or not
                # - this way there's the additional efficiency of
                #   handling things in bulk

                self._logger.info(
                    'Invoking run("rsync ...") data_path=%s tmpdir=%s...',
                    str(data_path.absolute()),
                    str(tmpdir),
                )
                self.update_task(ProcessState.PROCESSING, 'Resolving symlinks...')
                subprocess.run(
                    ["rsync", "-aL", str(data_path.absolute()) + "/", str(tmpdir)],
                    check=True,
                )

                self._logger.info(
                    'Invoking Popen("7z ...") tmpdir=%s...',
                    str(tmpdir),
                )
                self.update_task(ProcessState.PROCESSING, 'Compressing tarball...')
                compress_process = subprocess.Popen(
                    ["7z", "a", "-tzip", '-mmt=on', "-mx=4", tarball_path, "."],
                    cwd=str(tmpdir),
                )

                check_popen_progress(compress_process, poll_frequency)
                compress_process.wait()

                self._logger.info(
                    'Finished Popen("7z ...") tmpdir=%s returncode=%s',
                    str(tmpdir),
                    compress_process.returncode,
                )

        else:
            with open(tarball_path, "wb") as output_file:
                # due to the way gzip and pigz work, specifically, not
                # allowing to set the output file but use pipes instead, I
                # have to create 2 processes, one creates the zipped
                # stream and sends it to stdout, the other one catches it
                # and creates the file
                self.update_task(ProcessState.PROCESSING, 'Compressing tarball...')
                self._logger.info(
                    'Invoking Popen("tar ...") data_path=%s...',
                    data_path.absolute(),
                )
                tar_process = subprocess.Popen(
                    [
                        "tar",
                        "--dereference",
                        "--hard-dereference",
                        "-C",
                        data_path.absolute(),
                        "-cf",
                        "-",
                        ".",
                    ],
                    stdout=subprocess.PIPE,
                )
                self._logger.info(
                    'Invoking Popen("pigz ...") output=%s...',
                    tarball_path,
                )
                compress_process = subprocess.Popen(
                    ['pigz', '-4', "-c"],
                    stdin=tar_process.stdout,
                    stdout=output_file,
                )
                tar_process.stdout.close()  # type: ignore[union-attr]

                check_popen_progress(compress_process, poll_frequency)

                tar_process.wait()
                compress_process.wait()

                self._logger.info(
                    'Finished Popen("tar ...") data_path=%s name=%s returncode=%s',
                    str(data_path.absolute()),
                    data_path.name,
                    tar_process.returncode,
                )
                self._logger.info(
                    'Finished Popen("pigz ...") tarball_path=%s retruncode=%s',
                    tarball_path,
                    compress_process.returncode,
                )

        tarball_size_bytes: int = os.path.getsize(tarball_path)
        self._logger.info(
            "Created tarball %s (size=%s)",
            tarball_path,
            humanize.naturalsize(tarball_size_bytes, binary=True),
        )


def _write_text_pdf(pdf_path, message):
    """Write a minimal one-page PDF whose body is ``message``.

    Used as the README.pdf fallback when pandoc fails — built by hand so it
    has no dependency on pandoc/latex (which is exactly what's broken when
    we get here).
    """

    def _escape(text):
        return text.replace('\\', '\\\\').replace('(', '\\(').replace(')', '\\)')

    wrapped = []
    for raw in message.splitlines() or ['']:
        if not raw:
            wrapped.append('')
            continue
        while len(raw) > 80:
            wrapped.append(raw[:80])
            raw = raw[80:]
        wrapped.append(raw)

    text_ops = ['BT', '/F1 12 Tf', '50 750 Td', '14 TL']
    for line in wrapped:
        text_ops.append(f'({_escape(line)}) Tj T*')
    text_ops.append('ET')
    content = '\n'.join(text_ops).encode('latin-1', errors='replace')

    objects = [
        b'<< /Type /Catalog /Pages 2 0 R >>',
        b'<< /Type /Pages /Kids [3 0 R] /Count 1 >>',
        b'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] '
        b'/Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >>',
        b'<< /Length '
        + str(len(content)).encode('ascii')
        + b' >>\nstream\n'
        + content
        + b'\nendstream',
        b'<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>',
    ]

    out = bytearray(b'%PDF-1.4\n')
    offsets = []
    for i, obj in enumerate(objects, start=1):
        offsets.append(len(out))
        out += f'{i} 0 obj\n'.encode('ascii') + obj + b'\nendobj\n'

    xref_offset = len(out)
    out += f'xref\n0 {len(objects) + 1}\n'.encode('ascii')
    out += b'0000000000 65535 f \n'
    for off in offsets:
        out += f'{off:010d} 00000 n \n'.encode('ascii')
    out += (
        f'trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n'
        f'startxref\n{xref_offset}\n%%EOF\n'
    ).encode('ascii')

    with open(pdf_path, 'wb') as fh:
        fh.write(out)


def _is_mol_or_sdf(path):
    """Returns True if the file and path look like a MOL or SDF file.
    It does this by simply checking the file's extension.
    """
    return Path(path).suffix.lower() in ('.sdf', '.mol')


def _read_and_patch_molecule_name(path, molecule_name=None):
    """Patches the source file (expected to be a MOL or SDF file)
    by adding the molecule name (the basename of the file) and returning the
    file content as string. the assumption is that the source file is smll
    and can be read into memory.

    If a code/molecule is added we use that, otherwise we use the basename of
    the cleaned filename.

    Do not call this function for files that are not MOL or SD files.
    """
    logger.debug('Patching MOL/SDF "%s" molecule_name=%s', path, molecule_name)

    # The name will be set from file name
    # (without path prefix and the extension)
    # of the cleaned-up name.
    # e.g. the name of 'media/sdfs/Mpro-x3351_0A_rtEVbqf.sdf'
    # is 'Mpro-x3351_0A'.
    name = molecule_name or os.path.splitext(clean_filename(path))[0]

    # Now read the file, checking the first line
    # and setting it to the molecule name if it's blank.
    # We accumulate the file's content into 'content',
    # which we eventually return to the caller.
    content = ''
    with open(path, 'r', encoding='utf-8') as f_in:
        if first_line := f_in.readline().strip():
            content += first_line + '\n'
        else:
            content += name + '\n'
        # The rest of the file...
        for next_line in f_in:
            content += next_line

    # add sdf marker, the file read is mol but the combined file is sdf
    if Path(path).suffix.lower() != '.sdf':
        content += '$$$$\n\n'

    return content


def get_download_params(validated_data):
    """Extract download flags from serializer's validated data"""
    protein_params = {
        'pdb_info': validated_data['pdb_info'],
        'apo_file': validated_data['all_aligned_structures'],
        'bound_file': validated_data['all_aligned_structures'],
        'apo_solv_file': validated_data['all_aligned_structures'],
        'apo_desolv_file': validated_data['all_aligned_structures'],
        'ligand_pdb': validated_data['all_aligned_structures'],
        'ligand_sdf': validated_data['all_aligned_structures'],
        'ligand_smiles': validated_data['all_aligned_structures'],
        'cif_info': validated_data['cif_info'],
        'mtz_info': validated_data['mtz_info'],
        'map_info': validated_data['map_info'],
        'event_file': validated_data['event_file'],
        'sigmaa_file': validated_data['sigmaa_file'],
        'diff_file': validated_data['diff_file'],
    }

    other_params = {
        'sdf_info': validated_data['all_aligned_structures'],
        'single_sdf_file': validated_data['single_sdf_file'],
        'metadata_info': validated_data['metadata_info'],
        'smiles_info': validated_data['all_aligned_structures'],
        'trans_matrix_info': validated_data['trans_matrix_info'],
        'compound_sets': validated_data['compound_sets'],
        'soakdb_files': validated_data['soakdb_files'],
    }

    static_link = validated_data['static_link']

    return protein_params, other_params, static_link


def download_capacity_exceeded() -> bool:
    """Return True when too many download-build tasks are already in flight.

    The cap is derived dynamically from the celery workers' advertised
    concurrency: we sum each worker's pool max-concurrency, multiply by
    settings.MAX_DOWNLOAD_CONCURRENCY_PERCENT/100, and compare against the
    number of DownloadLinks rows that have a task_id but no file_url
    (i.e. still being built).
    """
    percent = settings.MAX_DOWNLOAD_CONCURRENCY_PERCENT
    if percent <= 0:
        return False

    try:
        worker_stats = celery_app.control.inspect().stats()
    except Exception as exc:
        logger.warning("Celery inspect.stats() failed (%s); skipping cap check", exc)
        return False

    if not worker_stats:
        logger.warning("No celery workers responded to stats(); skipping cap check")
        return False

    total_concurrency = 0
    for stats in worker_stats.values():
        try:
            total_concurrency += int(stats["pool"]["max-concurrency"])
        except (KeyError, TypeError, ValueError):
            continue

    if total_concurrency <= 0:
        logger.warning("Celery workers report zero concurrency; skipping cap check")
        return False

    cap = (total_concurrency * percent) // 100
    if cap <= 0:
        cap = 1

    in_progress = DownloadLinks.objects.filter(
        task_id__isnull=False,
        file_url__isnull=True,
    ).count()

    if in_progress >= cap:
        logger.warning(
            "Download capacity exceeded: in_progress=%d cap=%d "
            "(workers=%d total_concurrency=%d percent=%d)",
            in_progress,
            cap,
            len(worker_stats),
            total_concurrency,
            percent,
        )
        return True

    return False


def create_download(download_link_id: int, task, use_zip: bool = False):
    """Check/create a download zip file.

    This function is being ran inside a celery task, hence the object
    ids instead of objects themselves.

    Downloads are located in <MEDIA_ROOT>/downloads/ using a subdirectory
    using a UUID-4 value, with the file located in it, using the target title.
    For example: "/code/media/downloads/4c3afc69-bca9-4fb1-a76e-56c85a85899f/XX01ZVNS2B.zip".

    This function constructs the download file or returns a download form an exiting record.

    Returns:
        [file]: [URL to the file in the media directory]

    """

    # Retrieve the DOwnloadLinks record (and set its Task ID)
    download_link = DownloadLinks.objects.get(pk=download_link_id)
    download_link.task_id = str(task.request.id)
    download_link.save()

    # Simulate a worker crashing mid-build (task_id set, but file_url never written)
    # so the "lost task" housekeeping can be exercised. Ignored in production.
    if have_infection(INFECTION_STRUCTURE_DOWNLOAD_TASK):
        raise RuntimeError(f'Download Error! ({INFECTION_STRUCTURE_DOWNLOAD_TASK})')

    # error checking is not necessary because all these objects are
    # already resolved in the view and then passed through task
    site_observations = SiteObservation.objects.filter(
        pk__in=download_link.proteins.split(','),
    )

    # A custom logging adapter to refine the standard logger
    # by adding lightweight Task context to the messages we log
    _logger = TaskLoggerAdapter(
        logger,
        {
            'task': str(task.request.id),
            'marker': 'DOWNLOAD',
            'target': download_link.target,
            'tas': download_link.target.project,
            'username': download_link.user.username,
        },
    )

    # Always issue an INFO - we're one of the first functions to run
    # for this task. Here we provide the Task, Target, TAS and Username.
    # Other users of TaskLoggerAdapter can omit properties other than
    # the Task in order to shorten lines. The Task gives uas a key
    # correaltion value so that we can find matching lines.
    _logger.info('Handling download (target_id=%s)', download_link.target.pk)
    _logger.debug(
        'site_observation_ids=%s', site_observations.values_list('pk', flat=True)
    )

    task.update_state(
        state=ProcessState.PROCESSING,
        meta={
            "proposal_ref": download_link.target.project.title,
            "description": 'Start processing',
        },
    )

    if use_zip:
        filename = f'{download_link.target.title}.zip'
    else:
        filename = f'{download_link.target.title}.tar.gz'
    # The per-download directory uses the task_id (set above) so the path is
    # reproducible from the model alone — see DownloadLinks.get_file_url().
    file_url = os.path.join(
        settings.MEDIA_ROOT, 'downloads', download_link.task_id, filename
    )

    with TemporaryDirectory() as tempdir:
        downloader = DownloadStructures(
            task=task,
            tempdir=tempdir,
            target=download_link.target,
            use_zip=use_zip,
            target_access_string=download_link.target.project.title,
        )
        zip_contents = downloader.create_content_dict(
            site_obvs=site_observations,
            protein_params=download_link.protein_params,
            other_params=download_link.other_params,
        )
        downloader.create_tarball(
            zip_contents=zip_contents,
            file_url=file_url,
            original_search=download_link.original_search,
            site_observations=site_observations,
        )

    task.update_state(
        state=ProcessState.PROCESSING,
        meta={
            "proposal_ref": download_link.target.project.title,
            "description": 'File created',
        },
    )

    # We now have a download file
    # so record it and set the 'keep unitl' (expiry) time. file_url stores
    # the basename only; the absolute path is reconstructed by
    # DownloadLinks.get_file_url().
    download_link.file_url = filename
    download_link.keep_zip_until = datetime.now(timezone.utc) + KEEP_UNTIL_DURATION
    download_link.save()

    task.update_state(
        state=ProcessState.SUCCESS,
        meta={
            "proposal_ref": download_link.target.project.title,
            "description": file_url,
        },
    )
    _logger.info('- Handled new record (file_url=%s)', file_url)

    return file_url


def soft_erase_out_of_date_download_records():
    """Here we look for non-static records that are 'out of date' that have also
    not already been soft-deleted (i.e. where an 'expired_date' is not set)
    and we set the 'expired_date'.

    A record is considered out of date when its keep_zip_until is in the past.
    In-progress records (keep_zip_until still null) are handled separately by
    expire_lost_download_records(), which inspects the Celery task.

    Users should not use records that have an expired date.
    """
    now: datetime = datetime.now(timezone.utc)
    new_expired_records = (
        DownloadLinks.objects.filter(keep_zip_until__lt=now)
        .filter(expired_date__isnull=True)
        .filter(static_link=False)
    )

    if num_expired := new_expired_records.update(expired_date=now):
        msg = "1 record has" if num_expired == 1 else f"{num_expired} records have"
        logger.info('%s now expired', msg)


def expire_lost_download_records():
    """Find in-progress download records whose Celery task will never finish and
    expire them, recording a human-readable reason.

    An in-progress record is one that has not yet been given a keep_zip_until
    (that is only set when create_download() completes) and has not already been
    expired. Such a record can be left "stuck" forever if its task crashed or its
    worker was rebooted before it could write file_url / keep_zip_until.

    Each candidate is given TASK_START_GRACE_PERIOD to get going, then its task is
    inspected:

    - no task_id ........ the task was never launched -> lost
    - FAILURE / REVOKED . the task ended badly -> lost
    - PENDING ........... no worker picked it up, or the result was lost -> lost
    - still running ..... left alone, unless it has been running longer than
                          TASK_MAX_RUNTIME (a frozen result implies a lost worker)

    Users querying a lost task see the reason via TaskStatusView.
    """
    now: datetime = datetime.now(timezone.utc)
    start_cutoff = now - TASK_START_GRACE_PERIOD
    runtime_cutoff = now - TASK_MAX_RUNTIME

    in_progress_records = (
        DownloadLinks.objects.filter(keep_zip_until__isnull=True)
        .filter(expired_date__isnull=True)
        .filter(static_link=False)
    )

    num_lost = 0
    for record in in_progress_records:
        # Give the task a chance to start before we judge it.
        if record.create_date > start_cutoff:
            continue

        reason: str | None = None
        if not record.task_id:
            reason = "The download task was never launched"
        else:
            result = AsyncResult(str(record.task_id))
            state = result.state
            if state in ('FAILURE', 'REVOKED'):
                detail = result.info if isinstance(result.info, str) else state
                reason = f"The download task ended unexpectedly ({detail})"
            elif state == 'PENDING':
                reason = (
                    "The download task was lost"
                    " (no worker result - it may never have run)"
                )
            elif record.create_date < runtime_cutoff:
                runtime_m = settings.DOWNLOAD_TASK_MAX_RUNTIME_M
                reason = (
                    f"The download task exceeded the maximum runtime ({runtime_m} min)"
                    " and is presumed lost (e.g. the worker was restarted)"
                )

        if reason:
            record.expired_date = now
            record.expiry_reason = reason
            record.save(update_fields=['expired_date', 'expiry_reason'])
            logger.info(
                "Expired lost download (pk=%s task=%s): %s",
                record.pk,
                record.task_id,
                reason,
            )
            num_lost += 1

    if num_lost:
        msg = (
            "1 lost download has"
            if num_lost == 1
            else f"{num_lost} lost downloads have"
        )
        logger.info('%s now been expired', msg)


def hard_erase_out_of_date_download_records():
    """Physical zip files for non-static (dynamic) links
    are removed if they have been 'expired' for long enough.
    We can safely do this because the

    This is for security reasons and to conserve memory space.
    """
    now: datetime = datetime.now(timezone.utc)

    # Collect non-static records where there is an expired date but are not deleted
    # (where there is also a file-url)
    dead_dynamic_records = (
        DownloadLinks.objects.filter(expired_date__isnull=False)
        .filter(deleted__isnull=True)
        .filter(static_link=False)
        .filter(file_url__isnull=False)
    )

    num_deleted = 0
    num_pending = 0
    for dead_dynamic_record in dead_dynamic_records:
        # Skip any that have not been in an expired state for long enough...
        if dead_dynamic_record.expired_date + HARD_EXPIRY_GRACE_PERIOD >= now:
            num_pending += 1
            continue

        full_path = dead_dynamic_record.get_file_url()
        if full_path is None:
            # Already cleared (no file_url / task_id); nothing to remove.
            dead_dynamic_record.deleted = True
            dead_dynamic_record.save()
            continue
        dir_name = os.path.dirname(full_path)
        logger.debug('Deleting %s...', dir_name)
        if os.path.isdir(dir_name):
            try:
                shutil.rmtree(dir_name)
                logger.debug('Removed %s', dir_name)
            except Exception as ex:
                # Avoid potential race condition deleting the directory.
                # It's possible another call or stack Pod has deleted the directory
                # since we checked on line 1450. Issue a warning if the error
                # does not look like 'No such file or directory'...
                if 'No such file' not in str(ex):
                    logger.warning('Failed to remove %s (%s)', dir_name, ex)

        dead_dynamic_record.deleted = True
        dead_dynamic_record.save()
        num_deleted += 1

    if num_deleted:
        msg = "1 file has" if num_deleted == 1 else f"{num_deleted} files have"
        logger.info('%s now been deleted', msg)
    if num_pending:
        msg = "1 record has" if num_pending == 1 else f"{num_pending} records have"
        logger.info('%s expired (but not long enough for file deletion)', num_pending)


# TODO: issue with single_sdf file
def get_total_size(path: str) -> int:
    """Estimate data directory size (resolves symlinks)"""
    out = subprocess.check_output(["du", "-sbL", path], text=True)
    return int(out.split()[0])
