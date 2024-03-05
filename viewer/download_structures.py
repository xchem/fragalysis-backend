"""
download_structures.py

Methods for downloading a Target Zip file used by the download_structure API.
"""

import copy
import json
import logging
import os
import re
import shutil
import uuid
import zipfile
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from io import BytesIO, StringIO
from pathlib import Path
from typing import Any, Dict

import pandoc
from django.conf import settings
from django.db.models import Exists, OuterRef, Subquery

from viewer.models import (
    DownloadLinks,
    SiteObservation,
    SiteObservationTag,
    SiteObvsSiteObservationTag,
    TagCategory,
)
from viewer.utils import clean_filename

from .serializers import DownloadStructuresSerializer

logger = logging.getLogger(__name__)

# Length of time to keep records of dynamic links.
KEEP_UNTIL_DURATION = timedelta(minutes=90)

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

TAG_CATEGORIES = (
    'ConformerSites',
    'CanonSites',
    'CrystalformSites',
    'Quatassemblies',
    'Crystalforms',
)
CURATED_TAG_CATEGORIES = ('Series', 'Forum', 'Other')


class TagSubquery(Subquery):
    """Annotate SiteObservation with tag of given category"""

    def __init__(self, category):
        query = SiteObservationTag.objects.filter(
            pk=Subquery(
                SiteObvsSiteObservationTag.objects.filter(
                    site_observation=OuterRef('pk'),
                    site_obvs_tag__category=TagCategory.objects.get(
                        category=category,
                    ),
                ).values('site_obvs_tag')[:1]
            )
        ).values('tag')[0:1]
        super().__init__(query)


class CuratedTagSubquery(Exists):
    """Annotate SiteObservation with tag of given category"""

    def __init__(self, tag):
        query = SiteObvsSiteObservationTag.objects.filter(
            site_observation=OuterRef('pk'),
            site_obvs_tag=tag,
        )
        super().__init__(query)


@dataclass(frozen=True)
class ArchiveFile:
    path: str
    archive_path: str
    site_observation: SiteObservation | None = None


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
    },
    'molecules': {
        'sdf_files': {},
        'sdf_info': False,
        'single_sdf_file': False,
        'smiles_info': {},
    },
    'metadata_info': None,
    'trans_matrix_info': None,
}


_ERROR_FILE = 'errors.csv'

# unlike v1, metadata doesn't exist anymore, needs compiling
_METADATA_FILE = 'metadata.csv'


def _is_mol_or_sdf(path):
    """Returns True if the file and path look like a MOL or SDF file.
    It does this by simply checking the file's extension.
    """
    return Path(path).suffix.lower() in ('.sdf', '.mol')


def _add_empty_file(ziparchive, archive_path):
    """When file is missing, add an empty file to the archive.


    Used to send an explicit signal to the downloader that the file is
    missing.
    """
    logger.debug('+_add_empty_file: %s', archive_path)
    ziparchive.writestr(f'{archive_path}_FILE_NOT_IN_UPLOAD', BytesIO(b'').getvalue())


def _read_and_patch_molecule_name(path, molecule_name=None):
    """Patches the source file (expected to be a MOL or SDF file)
    by adding the molecule name (the basename of the file) and returning the
    file content as string. the assumption is that the source file is smll
    and can be read into memory.

    If a code/molecule is added we use that, otherwise we use the basename of
    the cleaned filename.

    Do not call this function for files that are not MOL or SD files.
    """
    assert _is_mol_or_sdf(path)

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

    return content


def _patch_molecule_name(site_observation):
    """Patch the MOL or SDF file with molecule name.

    Processes the content of ligand_mol attribute of the
    site_observation object. Returns the content as string.

    Alternative to _read_and_patch_molecule_name function above
    which operates on files. As ligand_mol is now stored as text,
    slightly different approach was necessary.

    """
    logger.debug('Patching MOL/SDF of "%s"', site_observation)

    # Now read the file, checking the first line
    # and setting it to the molecule name if it's blank.
    lines = site_observation.ligand_mol_file.split('\n')
    if not lines[0].strip():
        lines[0] = site_observation.long_code
    return '\n'.join(lines)


def _add_file_to_zip_aligned(ziparchive, code, archive_file):
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
    logger.debug('+_add_file_to_zip_aligned: %s, %s', code, archive_file)
    if not archive_file:
        # Odd - assume success
        logger.error('No filepath value')
        return True

    # calling str on archive_file.path because could be None
    filepath = str(Path(settings.MEDIA_ROOT).joinpath(str(archive_file.path)))
    if Path(filepath).is_file():
        if _is_mol_or_sdf(filepath):
            # It's a MOL or SD file.
            # Read and (potentially) adjust the file
            # and add to the archive as a string.
            content = _read_and_patch_molecule_name(filepath, molecule_name=code)
            ziparchive.writestr(archive_file.archive_path, content)
        else:
            # Copy the file without modification
            ziparchive.write(filepath, archive_file.archive_path)
        return True
    elif archive_file.site_observation:
        ziparchive.writestr(
            archive_file.archive_path,
            _patch_molecule_name(archive_file.site_observation),
        )
        return True
    else:
        logger.warning('filepath "%s" is not a file', filepath)
        _add_empty_file(ziparchive, archive_file.archive_path)

    return False


def _add_file_to_sdf(combined_sdf_file, archive_file):
    """Append the requested sdf file to the single sdf file provided.

    Args:
        combined_sdf: Handle of combined_sdf_file
        filepath: filepath from record

    Returns:
        [boolean]: [True of record added]
    """
    if not archive_file.path:
        # Odd - assume success
        logger.error('No filepath value')
        return True

    if archive_file.path and archive_file.path != 'None':
        with open(combined_sdf_file, 'a', encoding='utf-8') as f_out:
            patched_sdf_content = _patch_molecule_name(archive_file.site_observation)
            f_out.write(patched_sdf_content)
        return True
    else:
        logger.warning('filepath "%s" is not a file', archive_file.path)

    return False


def _protein_files_zip(zip_contents, ziparchive, error_file):
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
                if not _add_file_to_zip_aligned(ziparchive, prot, f):
                    error_file.write(f'{param},{prot},{f.archive_path}\n')
                    prot_errors += 1

    return prot_errors


def _molecule_files_zip(zip_contents, ziparchive, combined_sdf_file, error_file):
    """Write molecule (SD file) related data to the ZIP file
    Returns molecule errors
    """

    mol_errors = 0
    logger.info(
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
        ] is True and not _add_file_to_zip_aligned(
            ziparchive, prot.split(":")[0], archive_file
        ):
            error_file.write(f'sdf_info,{prot},{archive_file.path}\n')
            mol_errors += 1

        # Append sdf file on the Molecule record to the combined_sdf_file.
        if zip_contents['molecules'][
            'single_sdf_file'
        ] is True and not _add_file_to_sdf(combined_sdf_file, archive_file):
            error_file.write(f'single_sdf_file,{prot},{archive_file.path}\n')
            mol_errors += 1

    return mol_errors


def _smiles_files_zip(zip_contents, ziparchive, download_path):
    """Create and write the smiles file to the ZIP file"""
    smiles_filename = os.path.join(download_path, 'smiles.smi')
    logger.info('Creating SMILES file "%s"...', smiles_filename)

    num_smiles = 0
    with open(smiles_filename, 'w', encoding='utf-8') as smilesfile:
        for smi in zip_contents['molecules']['smiles_info']:
            logger.debug('Adding "%s"...', smi)
            smilesfile.write(f'{smi},')
            num_smiles += 1

    logger.info('Added %s SMILES', num_smiles)

    ziparchive.write(
        smiles_filename,
        os.path.join(_ZIP_FILEPATHS['smiles_info'], os.path.basename(smiles_filename)),
    )
    os.remove(smiles_filename)


def _trans_matrix_files_zip(ziparchive, target):
    """Add transformation matrices to archive.

    Note that this will always be the latest information - even for
    preserved searches.
    """
    logger.info('+ Processing trans matrix files')

    # grab the last set of files for this target
    experiment_upload = target.experimentupload_set.order_by('commit_datetime').last()

    trans_matrix_files = (
        experiment_upload.neighbourhood_transforms,
        experiment_upload.conformer_site_transforms,
        experiment_upload.reference_structure_transforms,
    )
    for tmf in trans_matrix_files:
        filepath = Path(settings.MEDIA_ROOT).joinpath(str(tmf))
        archive_path = os.path.join(
            _ZIP_FILEPATHS['trans_matrix_info'],
            Path(str(tmf)).name,
        )
        if filepath.is_file():
            ziparchive.write(
                filepath,
                archive_path,
            )
        else:
            logger.warning('File %s does not exist', Path(str(tmf)).name)
            _add_empty_file(ziparchive, archive_path)


def _metadate_file_zip(ziparchive, target):
    """Compile and add metadata file to archive."""
    logger.info('+ Processing metadata')

    annotations = {}
    values = ['code', 'longcode', 'cmpd__compound_code', 'smiles']
    header = ['Code', 'Long code', 'Compound code', 'Smiles']

    for category in TagCategory.objects.filter(category__in=TAG_CATEGORIES):
        tag = f'tag_{category.category.lower()}'
        values.append(tag)
        header.append(category.category)
        annotations[tag] = TagSubquery(category.category)

    pattern = re.compile(r'\W+')
    for tag in SiteObservationTag.objects.filter(
        category__in=TagCategory.objects.filter(category__in=CURATED_TAG_CATEGORIES),
        target=target,
    ):
        # for reasons unknown, mypy thinks tag is a string
        tagname = f'tag_{pattern.sub("_", tag.tag).strip().lower()}'  # type: ignore[attr-defined]
        values.append(tagname)
        header.append(f'[{tag.category}] {tag.tag}')  # type: ignore[attr-defined]
        annotations[tagname] = CuratedTagSubquery(tag)

    # fmt: off
    qs = SiteObservation.filter_manager.by_target(
        target=target,
    ).prefetch_related(
        'cmpd',
        'siteobservationtags',
    ).annotate(**annotations).values_list(*values)
    # fmt: on

    buff = StringIO()
    buff.write(','.join(header))
    buff.write('\n')
    for so_values in qs:
        buff.write(
            ','.join(
                [
                    str(k) if k else 'False' if isinstance(k, bool) else ''
                    for k in so_values
                ]
            )
        )
        buff.write('\n')

    ziparchive.writestr(_METADATA_FILE, buff.getvalue())
    logger.info('+ Processing metadata')


def _extra_files_zip(ziparchive, target):
    """If an extra info folder exists at the target root level, then
    copy the contents to the output file as is.
    Note that this will always be the latest information - even for
    preserved searches.
    """

    num_processed = 0
    num_extra_dir = 0
    for experiment_upload in target.experimentupload_set.order_by('commit_datetime'):
        extra_files = (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(experiment_upload.task_id)
        )

        # taking the latest upload for now
        # add unpacked zip directory
        extra_files = [d for d in list(extra_files.glob("*")) if d.is_dir()][0]

        # add upload_[d] dir
        extra_files = next(extra_files.glob("upload_*"))
        extra_files = extra_files.joinpath('extra_files')

        logger.debug('extra_files path 2: %s', extra_files)
        logger.info('Processing extra files (%s)...', extra_files)

        if extra_files.is_dir():
            num_extra_dir = num_extra_dir + 1
            for dirpath, _, files in os.walk(extra_files):
                for file in files:
                    filepath = os.path.join(dirpath, file)
                    logger.info('Adding extra file "%s"...', filepath)
                    ziparchive.write(
                        filepath,
                        os.path.join(
                            _ZIP_FILEPATHS[f'extra_files_{num_extra_dir}'], file
                        ),
                    )
                    num_processed += 1
        else:
            logger.info('Directory does not exist (%s)...', extra_files)

    if num_processed == 0:
        logger.info('No extra files found')
    else:
        logger.info('Processed %s extra files', num_processed)


def _yaml_files_zip(ziparchive, target, transforms_requested: bool = False) -> None:
    """Add all yaml files (except transforms) from upload to ziparchive"""

    for experiment_upload in target.experimentupload_set.order_by('commit_datetime'):
        yaml_paths = (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(experiment_upload.task_id)
        )

        transforms = [
            Path(f.name).name
            for f in (
                experiment_upload.neighbourhood_transforms,
                experiment_upload.neighbourhood_transforms,
                experiment_upload.neighbourhood_transforms,
            )
        ]
        # taking the latest upload for now
        # add unpacked zip directory
        yaml_paths = [d for d in list(yaml_paths.glob("*")) if d.is_dir()][0]

        # add upload_[d] dir
        yaml_paths = next(yaml_paths.glob("upload_*"))

        archive_path = Path('yaml_files').joinpath(yaml_paths.parts[-1])

        yaml_files = [
            f
            for f in list(yaml_paths.glob("*.yaml"))
            if f.is_file() and f.name not in transforms
        ]

        logger.info('Processing yaml files (%s)...', yaml_files)

        for file in yaml_files:
            logger.info('Adding yaml file "%s"...', file)
            if not transforms_requested and file.name == 'neighbourhoods.yaml':
                # don't add this file if transforms are not requested
                continue
            ziparchive.write(file, str(Path(archive_path).joinpath(file.name)))


def _document_file_zip(ziparchive, download_path, original_search, host):
    """Create the document file
    This consists of a template plus an added contents description.
    """
    # Don't need...
    del host

    logger.info('Creating documentation...')

    template_file = os.path.join("/code/doc_templates", "download_readme_template.md")
    readme_filepath = os.path.join(download_path, 'Readme.md')
    with open(readme_filepath, "a", encoding="utf-8") as readme:
        _build_readme(readme, original_search, template_file, ziparchive)

    # Convert markdown to pdf file
    pdf_filepath = os.path.join(download_path, 'Readme.pdf')
    doc = pandoc.read(open(readme_filepath, "r", encoding="utf-8").read())
    pandoc.write(doc, file=pdf_filepath, format='latex', options=["--columns=72"])

    ziparchive.write(pdf_filepath, os.path.join(_ZIP_FILEPATHS['readme'], 'README.pdf'))
    os.remove(readme_filepath)
    os.remove(pdf_filepath)


def _build_readme(readme, original_search, template_file, ziparchive):
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
        logger.warning('Could not find template file (%s)', template_file)

    # Files Included
    list_of_files = ziparchive.namelist()
    readme.write("\n## Files included\n")
    list_of_files.sort()
    for filename in list_of_files:
        readme.write(f'- {filename}' + '\n')


def _create_structures_zip(target, zip_contents, file_url, original_search, host):
    """Write a ZIP file containing data from an input dictionary."""

    logger.info('+ _create_structures_zip(%s)', target.title)
    logger.info('file_url="%s"', file_url)
    logger.info('single_sdf_file="%s"', zip_contents['molecules']['single_sdf_file'])
    logger.info('sdf_files=%s', zip_contents['molecules']['sdf_files'])

    logger.debug('zip_contents=%s', zip_contents)

    download_path = os.path.dirname(file_url)
    logger.info('Creating download path (%s)', download_path)
    os.makedirs(download_path, exist_ok=True)

    error_filename = os.path.join(download_path, _ERROR_FILE)
    error_file = open(error_filename, "w", encoding="utf-8")
    error_file.write("Param,Code,File not found when assembling download\n")
    errors = 0

    # If a single sdf file is also wanted then create file to
    # add sdf files to a file called {target}_combined.sdf.
    combined_sdf_file = None
    if zip_contents['molecules']['single_sdf_file'] is True:
        combined_sdf_file = os.path.join(download_path, f'{target.title}_combined.sdf')
        logger.info('combined_sdf_file=%s', combined_sdf_file)

    with zipfile.ZipFile(file_url, 'w', zipfile.ZIP_DEFLATED) as ziparchive:
        # Read through zip_contents to compile the file
        errors += _protein_files_zip(zip_contents, ziparchive, error_file)
        if errors > 0:
            logger.warning('After _protein_files_zip() errors=%s', errors)

        if zip_contents['molecules']['sdf_files']:
            errors_before = errors
            errors += _molecule_files_zip(
                zip_contents, ziparchive, combined_sdf_file, error_file
            )
            if errors > errors_before:
                logger.warning('After _molecule_files_zip() errors=%s', errors)

        # Add combined_sdf_file to the archive?
        if (
            zip_contents['molecules']['single_sdf_file'] is True
            and combined_sdf_file
            and os.path.isfile(combined_sdf_file)
        ):
            logger.info('Adding combined_sdf_file "%s"...', combined_sdf_file)
            ziparchive.write(
                combined_sdf_file,
                os.path.join(
                    _ZIP_FILEPATHS['single_sdf_file'],
                    os.path.basename(combined_sdf_file),
                ),
            )
            os.remove(combined_sdf_file)

        # If smiles info is required, then write one column for each molecule
        # to a smiles.smi file and then add to the archive.
        if zip_contents['molecules']['smiles_info']:
            _smiles_files_zip(zip_contents, ziparchive, download_path)

        # compile and add metadata.csv
        if zip_contents['metadata_info']:
            _metadate_file_zip(ziparchive, target)

        if zip_contents['trans_matrix_info']:
            _trans_matrix_files_zip(ziparchive, target)

        _extra_files_zip(ziparchive, target)

        _yaml_files_zip(
            ziparchive, target, transforms_requested=zip_contents['trans_matrix_info']
        )

        _document_file_zip(ziparchive, download_path, original_search, host)

        error_file.close()
        if errors > 0:
            logger.warning('errors=%s Adding %s to ziparchive', errors, _ERROR_FILE)
            ziparchive.write(error_filename, _ERROR_FILE)
        os.remove(error_filename)


def _protein_garbage_filter(proteins):
    """Garbage filter. It seems that Mpro has had some 'references_' added
    that are not being cleared up properly. Will look at this in future
    epic, but for now remove them from the download.

    Args:
        proteins

    Returns:
        [list]: [update protein list]
    """
    return proteins.exclude(code__startswith=r'references_')


def _create_structures_dict(site_obvs, protein_params, other_params):
    """Write a ZIP file containing data from an input dictionary

    Args:
        target
        proteins
        protein_params
        other_params

    Returns:
        [dict]: [dictionary containing the file contents]
    """
    logger.info('Processing %d SiteObservations', site_obvs.count())

    # Read through zip_params to compile the parameters
    zip_contents: Dict[str, Any] = copy.deepcopy(zip_template)
    for so in site_obvs:
        for param in protein_params:
            if protein_params[param] is True:
                if param in ['pdb_info', 'mtz_info', 'cif_info', 'map_info']:
                    # experiment object
                    model_attr = getattr(so.experiment, param)
                    logger.debug(
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
                        # here the model_attr is already stringified
                        try:
                            exp_path = re.search(r"x\d*", so.code).group(0)  # type: ignore[union-attr]
                        except AttributeError:
                            logger.error('Unexpected shortcodeformat: %s', so.code)
                            exp_path = so.code

                        apath = Path('crystallographic_files').joinpath(exp_path)
                        if model_attr and model_attr != 'None':
                            archive_path = str(
                                apath.joinpath(
                                    Path(f)
                                    .parts[-1]
                                    .replace(so.experiment.code, so.code)
                                )
                            )
                        else:
                            archive_path = str(apath.joinpath(param))
                        afile.append(ArchiveFile(path=f, archive_path=archive_path))

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
                    'diff_file',
                ]:
                    # siteobservation object

                    model_attr = getattr(so, param)
                    logger.debug(
                        'Adding param to zip: %s, value: %s', param, model_attr
                    )
                    apath = Path('aligned_files').joinpath(so.code)
                    if model_attr and model_attr != 'None':
                        archive_path = str(
                            apath.joinpath(
                                Path(model_attr.name)
                                .parts[-1]
                                .replace(so.longcode, so.code)
                            )
                        )
                    else:
                        archive_path = str(apath.joinpath(param))

                    afile = [
                        ArchiveFile(
                            path=model_attr.name,
                            archive_path=archive_path,
                        )
                    ]
                else:
                    logger.warning('Unexpected param: %s', param)
                    continue

                zip_contents['proteins'][param][so.code] = afile

    zip_contents['molecules']['single_sdf_file'] = other_params['single_sdf_file']
    zip_contents['molecules']['sdf_info'] = other_params['sdf_info']

    # sdf information is held as a file on the Molecule record.
    if other_params['sdf_info'] or other_params['single_sdf_file']:
        num_molecules_collected = 0
        num_missing_sd_files = 0
        for so in site_obvs:
            if so.ligand_mol_file:
                # There is an SD file (normal)
                # sdf info is now kept as text in db field
                archive_path = str(
                    Path('aligned_files').joinpath(so.code).joinpath(f'{so.code}.sdf')
                )
                # path is ignored when writing sdfs but mandatory field
                zip_contents['molecules']['sdf_files'].update(
                    {
                        ArchiveFile(
                            path=archive_path,
                            archive_path=archive_path,
                            site_observation=so,
                        ): so.code
                    }
                )
                num_molecules_collected += 1
            else:
                # No file value (odd).
                logger.warning(
                    "SiteObservation record's 'ligand_mol_file' isn't set (%s)", so
                )
                num_missing_sd_files += 1

        # Report (in the log) anomalies
        if num_molecules_collected == 0:
            logger.warning('No SD files collected')
        else:
            logger.info('%s SD files collected', num_molecules_collected)

        if site_obvs.count() != num_molecules_collected:
            logger.warning(
                'Expected %d files, got %d', site_obvs.count(), num_molecules_collected
            )

        if num_missing_sd_files > 0:
            logger.warning('%d missing files', num_missing_sd_files)

    # The smiles at molecule level may not be unique.
    if other_params['smiles_info'] is True:
        for molecule in site_obvs:
            zip_contents['molecules']['smiles_info'].update({molecule.smiles: None})

    zip_contents['metadata_info'] = other_params['metadata_info']

    # Add the trans matrix files
    zip_contents['trans_matrix_info'] = other_params['trans_matrix_info']

    return zip_contents


def get_download_params(request):
    """Check whether structures have been previously downloaded

    Args:
        request

    Returns:
        protein_params, other_params
    """

    serializer = DownloadStructuresSerializer(data=request.data)
    valid = serializer.is_valid()
    logger.debug('serializer validated data: %s, %s', valid, serializer.validated_data)
    if not valid:
        logger.error('serializer errors: %s', serializer.errors)

    protein_params = {
        'pdb_info': serializer.validated_data['pdb_info'],
        'apo_file': serializer.validated_data['all_aligned_structures'],
        'bound_file': serializer.validated_data['all_aligned_structures'],
        'apo_solv_file': serializer.validated_data['all_aligned_structures'],
        'apo_desolv_file': serializer.validated_data['all_aligned_structures'],
        'ligand_pdb': serializer.validated_data['all_aligned_structures'],
        'cif_info': serializer.validated_data['cif_info'],
        'mtz_info': serializer.validated_data['mtz_info'],
        'map_info': serializer.validated_data['map_info'],
        'event_file': serializer.validated_data['event_file'],
        'sigmaa_file': serializer.validated_data['sigmaa_file'],
        'diff_file': serializer.validated_data['diff_file'],
    }

    other_params = {
        'sdf_info': serializer.validated_data['all_aligned_structures'],
        'single_sdf_file': serializer.validated_data['single_sdf_file'],
        'metadata_info': serializer.validated_data['metadata_info'],
        'smiles_info': serializer.validated_data['all_aligned_structures'],
        'trans_matrix_info': serializer.validated_data['trans_matrix_info'],
    }

    static_link = serializer.validated_data['static_link']

    return protein_params, other_params, static_link


def create_or_return_download_link(request, target, site_observations):
    """Check/create a download zip file.

    Downloads are located in <MEDIA_ROOT>/downloads/ using a subdirectory
    using a UUID-4 value, with the file located in it, using the target title.
    For example: "/code/media/downloads/4c3afc69-bca9-4fb1-a76e-56c85a85899f/XX01ZVNS2B.zip".

    This function constructs the download file or returns a download form an exiting record.

    Returns:
        [file]: [URL to the file in the media directory]
    """
    logger.info('+ Handling download for Target "%s"', target.title)

    # Log the provided SiteObservations
    num_given_site_obs = site_observations.count()
    site_ob_repr = "".join(
        # & syntax copies the queryset without evaluating it. this way
        # I have an unsliced queryset for later for protein_garbage
        # filter method (this is what gave sliced queryset error)
        "%r " % site_ob
        for site_ob in site_observations & site_observations
    )
    logger.debug(
        'Given %s SiteObservation records: %r', num_given_site_obs, site_ob_repr
    )

    protein_params, other_params, static_link = get_download_params(request)
    logger.debug('proteins_params: %s', protein_params)
    logger.debug('other_params: %s', other_params)
    logger.debug('static_link: %s', static_link)

    # Remove 'references_' from protein list if present.
    site_observations = _protein_garbage_filter(site_observations)
    if num_given_site_obs > site_observations.count():
        logger.warning(
            'Removed %d "references_" proteins from download',
            num_given_site_obs - site_observations.count(),
        )

    # Save the list of protein codes - this is the ispybsafe set for this user.
    proteins_list = list(site_observations.values_list('code', flat=True))
    logger.debug('proteins_list: %s', proteins_list)

    # Remove the token so the original search can be stored
    original_search = copy.deepcopy(request.data)
    if 'csrfmiddlewaretoken' in original_search:
        del original_search['csrfmiddlewaretoken']

    if existing_link := DownloadLinks.objects.filter(
        target_id=target.id,
        proteins=proteins_list,
        protein_params=protein_params,
        other_params=other_params,
    ).first():
        # Dynamic to static?
        # Static link records are never removed.
        if static_link and not existing_link.static_link:
            logger.info(
                'Converting dynamic link to static link (%s)', existing_link.file_url
            )
            existing_link.static_link = True
            existing_link.save()
        # Now return the file...
        file_url = existing_link.file_url
        assert os.path.isfile(file_url)
        logger.info('- Handled existing download (file_url=%s)', file_url)
        return file_url

    # No existing Download record - create one,
    # which requires construction of the file prior to creating the record.
    # A record indicates the file is present. It is removed
    # when "out of date".
    host = request.get_host()
    filename = f'{target.title}.zip'
    file_url = os.path.join(
        settings.MEDIA_ROOT, 'downloads', str(uuid.uuid4()), filename
    )
    logger.info('Creating new download (file_url=%s)...', file_url)

    zip_contents = _create_structures_dict(
        site_observations, protein_params, other_params
    )
    _create_structures_zip(target, zip_contents, file_url, original_search, host)

    download_link = DownloadLinks()
    # Note: 'zip_file' and 'zip_contents' record properties are no longer used.
    download_link.file_url = file_url
    download_link.user = request.user if request.user.is_authenticated else None
    download_link.target = target
    download_link.proteins = proteins_list
    download_link.protein_params = protein_params
    download_link.other_params = other_params
    download_link.static_link = static_link
    download_link.create_date = datetime.now(timezone.utc)
    download_link.original_search = original_search
    # We've just created the file, so the download is valid now...
    # Dynamic files are typically removed on the next download request
    # that occurs after the KEEP_UNTIL_DURATION.
    download_link.keep_zip_until = download_link.create_date + KEEP_UNTIL_DURATION
    download_link.save()

    logger.info('- Handled new record (file_url=%s)', file_url)
    return file_url


def erase_out_of_date_download_records():
    """Physical zip files and DownloadLink records for non-static (dynamic) links
    are removed after 1 hour (typically during a POST call to create a new download).

    This is for security reasons and to conserve memory space. Only if the file can
    be deleted do we delete the download record. So, if there are any problems
    with the file-system the model should continue to reflect the current state
    of the world.
    """
    num_removed = 0
    out_of_date_dynamic_records = DownloadLinks.objects.filter(
        keep_zip_until__lt=datetime.now(timezone.utc)
    ).filter(static_link=False)
    for out_of_date_dynamic_record in out_of_date_dynamic_records:
        file_url = out_of_date_dynamic_record.file_url
        logger.info(
            '+ Attempting to remove download link record (file_url=%s)...', file_url
        )

        dir_name = os.path.dirname(file_url)
        if os.path.isdir(dir_name):
            logger.debug('Removing file_url directory (%s)...', dir_name)
            shutil.rmtree(dir_name, ignore_errors=True)
            logger.debug('Removed (%s)', dir_name)

        # Does the file exist now?
        # Hopefully not - but cater for 'cosmic-ray-effect' and
        # only delete the originating record if the file has been removed.
        if os.path.isdir(dir_name):
            logger.warning(
                'Failed removal of file_url directory (%s), leaving record',
                dir_name,
            )
        else:
            logger.info(
                'Removed file_url directory (%s), removing DownloadLinks record...',
                dir_name,
            )
            out_of_date_dynamic_record.delete()
            num_removed += 1

    logger.info('Erased %d', num_removed)
