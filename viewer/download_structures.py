"""
download_structures.py

Methods for downloading a Target Zip file used by the download_structure API.
"""

import copy
import datetime
import json
import logging
import os
import shutil
import uuid
import zipfile
from pathlib import Path
from typing import Any, Dict

import pandoc
from django.conf import settings

from viewer.models import DownloadLinks
from viewer.utils import clean_filename

logger = logging.getLogger(__name__)

# Filepaths mapping for writing associated files to the zip archive.
# Note that if this is set to 'aligned' then the files will be placed in
# the protein code subdirectory of the aligned directory
# (as for the target upload).
_ZIP_FILEPATHS = {
    'apo_file': ('aligned'),
    'bound_file': ('aligned'),
    'cif_info': ('aligned'),
    'mtz_info': ('aligned'),
    # 'map_info': ('aligned'),
    'sigmaa_file': ('aligned'),
    'diff_file': ('aligned'),
    'event_file': ('aligned'),
    'sdf_info': ('aligned'),
    'single_sdf_file': (''),
    'metadata_info': (''),
    'smiles_info': (''),
    'extra_files': ('extra_files'),
    'readme': (''),
}

# Dictionary containing all references needed to create the zip file
# NB you may need to add a version number to this at some point...
zip_template = {
    'proteins': {
        'apo_file': {},  # from experiment
        'bound_file': {},  # x
        'cif_info': {},  # from experiment
        'mtz_info': {},  # from experiment
        'event_file': {},  # x
        'diff_file': {},  # renamed from diff_file and sigmaa_file
        'sigmaa_file': {},
    },
    'molecules': {
        'sdf_files': {},
        'sdf_info': False,
        'single_sdf_file': False,
        'smiles_info': {},
    },
    'metadata_info': None,
}


# A directory, relative to the media directory,
# where missing SD files are written.
# The SD files are constructed from the molecule 'sdf_info' field
# (essentially MOL-file text) when the 'sdf_file' field is blank.
_MISSING_SDF_DIRECTORY = 'missing-sdfs'
_MISSING_SDF_PATH = os.path.join(settings.MEDIA_ROOT, _MISSING_SDF_DIRECTORY)

_ERROR_FILE = 'errors.csv'


def _replace_missing_sdf(molecule, code):
    """Creates a file in the 'missing SDFs' directory, using the protein code
    provided. The file is constructed using the molecule's sdf_info field, skipping the
    action if the file exists. The media-relative path of the written file is returned
    (if it was written).

    Files, once written, are left and are not removed (or replaced).
    The directory serves an archive of missing SD files.

    This was added for FE/915 to generate SD files for those that are missing
    from the upload directory.
    """
    if not os.path.isdir(_MISSING_SDF_PATH):
        os.mkdir(_MISSING_SDF_PATH)

    # We shouldn't be called if molecule['sdf_info'] is blank.
    # but check anyway.
    sdf_info = molecule.ligand_mol_file
    if not sdf_info:
        return None
    sdf_lines = sdf_info.splitlines(True)[1:]
    if not sdf_lines:
        return None
    # Make sure last line ends with a new-line
    if not sdf_lines[-1].endswith('\n'):
        sdf_lines[-1] += '\n'

    # media-relative path to missing file...
    missing_file = os.path.join(_MISSING_SDF_DIRECTORY, f'{code}.sdf')
    # absolute path to missing file...
    missing_path = os.path.join(settings.MEDIA_ROOT, missing_file)
    # create the file if it doesn't exist...
    if not os.path.isfile(missing_path):
        # No file - create one.
        with open(missing_path, 'w', encoding='utf-8') as sd_file:
            # First line is the protein code, i.e. "PGN_RS02895PGA-x0346_0B"
            sd_file.write(f'{code}\n')
            # Now write the lines from the molecule sdf_info record
            sd_file.writelines(sdf_lines)
            # And append file terminator...
            sd_file.write('$$$$\n')

    # Returns the media-relative path to the file in the missing file directory
    return missing_file


def _add_file_to_zip(ziparchive, param, filepath):
    """Add the requested file to the zip archive.

    Args:
        ziparchive: Handle of zip archive
        param: parameter of filelist
        filepath: filepath from record

    Returns:
        [boolean]: [True of record added to error file]
    """
    media_root = settings.MEDIA_ROOT
    if not filepath:
        return False

    fullpath = os.path.join(media_root, filepath)

    if os.path.isfile(fullpath):
        cleaned_filename = clean_filename(filepath)
        ziparchive.write(
            fullpath, os.path.join(_ZIP_FILEPATHS[param], cleaned_filename)
        )
        return False

    return True


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
    assert _is_mol_or_sdf(path)

    logger.debug('Patching MOL/SDF "%s" molecule_name=%s', path, molecule_name)

    name = molecule_name
    if not name:
        # No molecule name provided.
        # The name will be set from file name
        # (without path prefix and the extension)
        # of the cleaned-up name.
        # e.g. the name of 'media/sdfs/Mpro-x3351_0A_rtEVbqf.sdf'
        # is 'Mpro-x3351_0A'.
        name = os.path.splitext(clean_filename(path))[0]

    # Now read the file, checking the first line
    # and setting it to the molecule name if it's blank.
    # We accumulate the file's content into 'content',
    # which we eventually return to the caller.
    content = ''
    with open(path, 'r', encoding='utf-8') as f_in:
        # First line (stripped)
        first_line = f_in.readline().strip()
        if first_line:
            content += first_line + '\n'
        else:
            content += name + '\n'
        # The rest of the file...
        for next_line in f_in:
            content += next_line

    return content


def _add_file_to_zip_aligned(ziparchive, code, filepath):
    """Add the requested file to the zip archive.

    If the file is an SDF or MOL we insert the name of the molecule
    if it is not set  (the basename of the file).

    Args:
        ziparchive: Handle of zip archive
        code: original protein code stripped of any alternate name.
        filepath: filepath from record

    Returns:
        [boolean]: [True of record added to error file]
    """
    if not filepath:
        return False

    logger.debug('+ _add_file_to_zip_aligned, filepath: %s', filepath)

    # apparently incoming filepath can be both str and FieldFile
    try:
        filepath = filepath.path
    except AttributeError:
        filepath = str(Path(settings.MEDIA_ROOT).joinpath(filepath))

    if Path(filepath).is_file():
        archive_path = str(Path(*Path(filepath).parts[7:]))
        if _is_mol_or_sdf(filepath):
            # It's a MOL or SD file.
            # Read and (potentially) adjust the file
            # and add to the archive as a string.
            content = _read_and_patch_molecule_name(filepath, molecule_name=code)
            ziparchive.writestr(archive_path, content)
        else:
            # Copy the file without modification
            ziparchive.write(filepath, archive_path)
        return False

    return True


def _add_file_to_sdf(combined_sdf_file, filepath):
    """Append the requested sdf file to the single sdf file provided.

    Args:
        combined_sdf: Handle of combined_sdf_file
        filepath: filepath from record

    Returns:
        [boolean]: [True of record added to error file]
    """
    media_root = settings.MEDIA_ROOT

    if not filepath:
        return False

    fullpath = os.path.join(media_root, filepath)

    if os.path.isfile(fullpath):
        with open(combined_sdf_file, 'a', encoding='utf-8') as f_out:
            patched_sdf_content = _read_and_patch_molecule_name(fullpath)
            f_out.write(patched_sdf_content)
        return False

    return True


def _protein_files_zip(zip_contents, ziparchive, error_file):
    """Write all protein related data to the ZIP file
    Returns protein errors
    """
    prot_errors = 0
    for param, files in zip_contents['proteins'].items():
        if not files:
            continue

        for prot, prot_file in files.items():
            logger.debug('%s: %s', prot, prot_file)
            if _add_file_to_zip_aligned(ziparchive, prot.split(":")[0], prot_file):
                error_file.write('{},{},{}\n'.format(param, prot, prot_file))
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
    for file, prot in zip_contents['molecules']['sdf_files'].items():
        # Do not try and process any missing SD files.
        if not file:
            error_file.write('{},{},{}\n'.format('sdf_files', prot, 'missing'))
            mol_errors += 1
            continue

        if zip_contents['molecules']['sdf_info'] is True:
            # Add sdf file on the Molecule record to the archive folder.
            if _add_file_to_zip_aligned(ziparchive, prot.split(":")[0], file):
                error_file.write('{},{},{}\n'.format('sdf_info', prot, file))
                mol_errors += 1

        # Append sdf file on the Molecule record to the combined_sdf_file.
        if zip_contents['molecules']['single_sdf_file'] is True:
            if _add_file_to_sdf(combined_sdf_file, file):
                error_file.write('{},{},{}\n'.format('single_sdf_file', prot, file))
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
            smilesfile.write(smi + ',')
            num_smiles += 1

    logger.info('Added %s SMILES', num_smiles)

    ziparchive.write(
        smiles_filename,
        os.path.join(_ZIP_FILEPATHS['smiles_info'], os.path.basename(smiles_filename)),
    )
    os.remove(smiles_filename)


def _extra_files_zip(ziparchive, target):
    """If an extra info folder exists at the target root level, then
    copy the contents to the output file as is.
    Note that this will always be the latest information - even for
    preserved searches.
    """

    num_processed = 0
    extra_files = os.path.join(
        settings.MEDIA_ROOT, 'targets', target.title, 'extra_files'
    )
    logger.info('Processing extra files (%s)...', extra_files)

    if os.path.isdir(extra_files):
        for dirpath, _, files in os.walk(extra_files):
            for file in files:
                filepath = os.path.join(dirpath, file)
                logger.info('Adding extra file "%s"...', filepath)
                ziparchive.write(
                    filepath, os.path.join(_ZIP_FILEPATHS['extra_files'], file)
                )
                num_processed += 1
    else:
        logger.info('Directory does not exist (%s)...', extra_files)

    if num_processed == 0:
        logger.info('No extra files found')
    else:
        logger.info('Processed %s extra files', num_processed)


# def _get_external_download_url(download_path, host):
#     """Returns the external download URL from the internal url for
#     the documentation.
#     This a bit messy but requirements change and this should be replaced
#     by data from the frontend in a future issue.
#     """
#
#     download_base = os.path.join(settings.MEDIA_ROOT, 'downloads')
#     download_uuid = download_path.replace(download_base, "")
#     external_path = os.path.join(
#         settings.SECURE_PROXY_SSL_HEADER[1] + '://' + host,
#         'viewer/react/download/tag')
#     external_path = external_path+download_uuid
#
#     return external_path


def _document_file_zip(ziparchive, download_path, original_search, host):
    """Create the document file
    This consists of a template plus an added contents description.
    """
    # Don't need...
    del host

    logger.info('Creating documentation...')

    template_file = os.path.join("/code/doc_templates", "download_readme_template.md")

    readme_filepath = os.path.join(download_path, 'Readme.md')
    pdf_filepath = os.path.join(download_path, 'Readme.pdf')

    with open(readme_filepath, "a", encoding="utf-8") as readme:
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
        readme.write("```" + json.dumps(original_search) + "```\n\n")

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
            readme.write('- ' + filename + '\n')

    # Convert markdown to pdf file
    doc = pandoc.read(open(readme_filepath, "r", encoding="utf-8").read())
    pandoc.write(doc, file=pdf_filepath, format='latex', options=["--columns=72"])

    ziparchive.write(pdf_filepath, os.path.join(_ZIP_FILEPATHS['readme'], 'README.pdf'))
    os.remove(readme_filepath)
    os.remove(pdf_filepath)


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
    error_file.write("Param,Code,Invalid file reference\n")
    errors = 0

    # If a single sdf file is also wanted then create file to
    # add sdf files to a file called {target}_combined.sdf.
    combined_sdf_file = None
    if zip_contents['molecules']['single_sdf_file'] is True:
        combined_sdf_file = os.path.join(
            download_path, '{}_combined.sdf'.format(target.title)
        )
        logger.info('combined_sdf_file=%s', combined_sdf_file)

    with zipfile.ZipFile(file_url, 'w') as ziparchive:
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

        # Add the metadata file from the target
        if zip_contents['metadata_info']:
            if _add_file_to_zip(
                ziparchive, 'metadata_info', zip_contents['metadata_info']
            ):
                error_file.write(
                    '{},{},{}\n'.format(
                        'metadata_info', target, zip_contents['metadata_info']
                    )
                )
                errors += 1
                logger.warning('After _add_file_to_zip() errors=%s', errors)

        _extra_files_zip(ziparchive, target)

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


def _create_structures_dict(target, site_obvs, protein_params, other_params):
    """Write a ZIP file containing data from an input dictionary

    Args:
        target
        proteins
        protein_params
        other_params

    Returns:
        [dict]: [dictionary containing the file contents]
    """
    logger.info('+ _create_structures_dict')

    zip_contents: Dict[str, Any] = copy.deepcopy(zip_template)

    logger.info(
        'Expecting %d molecules from %d proteins', site_obvs.count(), site_obvs.count()
    )

    # Read through zip_params to compile the parameters
    for so in site_obvs:
        for param in protein_params:
            if protein_params[param] is True:
                try:
                    # getting the param from experiment. more data are
                    # coming from there, that's why this is in try
                    # block
                    # getattr retrieves FieldFile object, hance the .name
                    zip_contents['proteins'][param][so.code] = getattr(
                        so.experiment, param
                    ).name
                except AttributeError:
                    # on the off chance that the data are in site_observation model
                    zip_contents['proteins'][param][so.code] = getattr(so, param).name

    if other_params['single_sdf_file'] is True:
        zip_contents['molecules']['single_sdf_file'] = True

    if other_params['sdf_info'] is True:
        zip_contents['molecules']['sdf_info'] = True

    # sdf information is held as a file on the Molecule record.
    num_molecules_collected = 0
    if other_params['sdf_info'] or other_params['single_sdf_file']:
        num_missing_sd_files = 0
        for so in site_obvs:
            rel_sd_file = None
            if so.ligand_mol_file:
                # There is an SD file (normal)
                # sdf info is now kept as text in db field
                rel_sd_file = _replace_missing_sdf(so, so.code)
            else:
                # No file value (odd).
                logger.warning(
                    "SiteObservation record's 'ligand_mol_file' isn't set (protein.code=%s).",
                    so.code,
                )
                logger.error(
                    "Molecule record's 'sdf_info' isn't set (protein.code=%s).", so.code
                )
                num_missing_sd_files += 1

            if rel_sd_file:
                logger.debug('rel_sd_file=%s protein.code=%s', rel_sd_file, so.code)
                zip_contents['molecules']['sdf_files'].update({rel_sd_file: so.code})
                num_molecules_collected += 1

        # Report (in the log) anomalies
        if num_molecules_collected == 0:
            logger.warning('No SD files collected')
        else:
            logger.info('%s SD files collected', num_molecules_collected)

        if site_obvs.count() != num_molecules_collected:
            logger.error('Expected %d files', site_obvs.count())

        if num_missing_sd_files > 0:
            logger.error('%d missing files', num_missing_sd_files)

    # The smiles at molecule level may not be unique.
    if other_params['smiles_info'] is True:
        for molecule in site_obvs:
            zip_contents['molecules']['smiles_info'].update({molecule.smiles: None})

    # Add the metadata file from the target
    if other_params['metadata_info'] is True:
        zip_contents['metadata_info'] = target.metadata.name

    return zip_contents


def get_keep_until():
    """Return the keep_until time"""
    return datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=1)


def get_download_params(request):
    """Check whether structures have been previously downloaded

    Args:
        request

    Returns:
        protein_params, other_params
    """
    protein_param_flags = [
        'apo_file',
        'bound_file',
        'cif_info',
        'mtz_info',
        'event_file',
        'sigmaa_file',
        'diff_file',
    ]

    other_param_flags = ['sdf_info', 'single_sdf_file', 'metadata_info', 'smiles_info']

    # protein_params = {'pdb_info': request.data['pdb_info'],
    #               'bound_info': request.data['bound_info'],
    #               'cif_info': request.data['cif_info'],
    #               'mtz_info': request.data['mtz_info'],
    #               'diff_info': request.data['diff_info'],
    #               'event_info': request.data['event_info'],
    #               'sigmaa_info': request.data['sigmaa_info'],
    #               'trans_matrix_info':
    #                   request.data['trans_matrix_info']}
    protein_params = {}
    for param in protein_param_flags:
        protein_params[param] = False
        if param in request.data:
            if request.data[param] == True or request.data[param] == 'true':
                protein_params[param] = True

    # other_params = {'sdf_info': request.data['sdf_info'],
    #                 'single_sdf_file': request.data['single_sdf_file'],
    #                 'metadata_info': request.data['metadata_info'],
    #                 'smiles_info': request.data['smiles_info']}
    other_params = {}
    for param in other_param_flags:
        other_params[param] = False
        if param in request.data:
            if request.data[param] == True or request.data[param] == 'true':
                other_params[param] = True

    static_link = False
    if 'static_link' in request.data:
        if request.data['static_link'] is True or request.data['static_link'] == 'true':
            static_link = True

    return protein_params, other_params, static_link


def check_download_links(request, target, site_observations):
    """Check/create the download zip file for dynamic links

    Args:
        request
        target
        proteins

    Returns:
        [file]: [URL to the file in the media directory]
    """

    logger.info('+ check_download_links(%s)', target.title)
    host = request.get_host()

    protein_params, other_params, static_link = get_download_params(request)
    logger.debug('proteins_params: %s', protein_params)
    logger.debug('other_params: %s', other_params)
    logger.debug('static_link: %s', static_link)

    # Remove 'references_' from protein list if present.
    num_given_proteins = site_observations.count()
    site_observations = _protein_garbage_filter(site_observations)
    num_removed = num_given_proteins - site_observations.count()
    if num_removed:
        logger.warning('Removed %d "references_" proteins from download', num_removed)

    # Save the list of protein codes - this is the ispybsafe set for
    # this user.
    # proteins_list = [p.code for p in site_observations]
    proteins_list = list(site_observations.values_list('code', flat=True))
    logger.debug('proteins_list: %s', proteins_list)

    # Remove the token so the original search can be stored
    original_search = copy.deepcopy(request.data)
    if 'csrfmiddlewaretoken' in original_search:
        del original_search['csrfmiddlewaretoken']

    # For dynamic files, if the target, proteins, protein_params and other_params
    # are in an existing record then this is a repeat search .
    # If the zip is there, then we can just return the file_url.
    # If the record is there but the zip file isn't, then regenerate the zip file.
    # from the search, to contain the latest information.
    # If no record is not there at all, then this is a new link.

    existing_link = DownloadLinks.objects.filter(
        target_id=target.id,
        proteins=proteins_list,
        protein_params=protein_params,
        other_params=other_params,
        static_link=False,
    )

    if existing_link.exists():
        if (
            existing_link[0].zip_file
            and os.path.isfile(existing_link[0].file_url)
            and not static_link
        ):
            logger.info('Existing link (%s)', existing_link[0].file_url)
            return existing_link[0].file_url, True

        if os.path.isfile(existing_link[0].file_url) and not static_link:
            logger.info('Existing link in progress (%s)', existing_link[0].file_url)
            # Repeat call, file is currently being created by another process.
            return existing_link[0].file_url, False

        # Recreate file.
        #
        # This step can result in references to missing SD Files (see FE issue 907).
        # If so the missing file will have a file reference of 'MISSING/<filename>'
        # in the corresponding ['molecules']['sdf_files'] entry.
        zip_contents = _create_structures_dict(
            target, site_observations, protein_params, other_params
        )
        _create_structures_zip(
            target,
            zip_contents,
            existing_link[0].file_url,
            existing_link[0].original_search,
            host,
        )
        existing_link[0].keep_zip_until = get_keep_until()
        existing_link[0].zip_file = True
        if static_link:
            # Changing from a dynamic to a static link
            existing_link[0].static_link = static_link
            existing_link[0].zip_contents = zip_contents
        existing_link[0].save()
        logger.info('Link saved (%s)', existing_link[0].file_url)

        return existing_link[0].file_url, True

    # Create a new link for this user
    # /code/media/downloads/<download_uuid>
    filename = target.title + '.zip'
    media_root = settings.MEDIA_ROOT
    file_url = os.path.join(media_root, 'downloads', str(uuid.uuid4()), filename)

    zip_contents = _create_structures_dict(
        target, site_observations, protein_params, other_params
    )
    _create_structures_zip(target, zip_contents, file_url, original_search, host)

    download_link = DownloadLinks()
    download_link.file_url = file_url
    if request.user.is_authenticated:
        download_link.user = request.user
    else:
        download_link.user = None
    download_link.target = target
    download_link.proteins = proteins_list
    download_link.protein_params = protein_params
    download_link.other_params = other_params
    download_link.static_link = static_link
    if static_link:
        download_link.zip_contents = zip_contents
    download_link.create_date = datetime.datetime.now(datetime.timezone.utc)
    download_link.keep_zip_until = get_keep_until()
    download_link.zip_file = True
    download_link.original_search = original_search
    download_link.save()

    return file_url, True


def recreate_static_file(existing_link, host):
    """Recreate static file for existing link"""
    target = existing_link.target

    _create_structures_zip(
        target,
        existing_link.zip_contents,
        existing_link.file_url,
        existing_link.original_search,
        host,
    )
    existing_link.keep_zip_until = get_keep_until()
    existing_link.zip_file = True
    existing_link.save()


def maintain_download_links():
    """Maintain zip files

    Physical zip files are removed after 1 hour (or when the next POST call
    is made) for security reasons and to conserve memory space.

    """

    old_links = DownloadLinks.objects.filter(
        keep_zip_until__lt=datetime.datetime.now(datetime.timezone.utc)
    ).filter(zip_file=True)
    for old_link in old_links:
        old_link.zip_file = False
        old_link.save()
        shutil.rmtree(os.path.dirname(old_link.file_url), ignore_errors=True)
