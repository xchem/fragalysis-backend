"""
download_structures.py

Methods for downloading a Target Zip file used by the download_structure API.
"""

import os
import zipfile
import datetime
import uuid
import shutil
import fnmatch
import logging
import copy

from django.conf import settings

from viewer.models import (
    Molecule,
    Protein,
    DownloadLinks
)

logger = logging.getLogger(__name__)

# Filepaths mapping for writing associated files to the zip archive.
# Note that if this is set to 'aligned' then the files will be placed in
# the protein code subdirectory of the aligned directory
# (as for the target upload).
_ZIP_FILEPATHS = {
    'pdb_info': ('aligned'),
    'bound_info': ('aligned'),
    'cif_info': ('aligned'),
    'mtz_info': ('aligned'),
    'map_info': ('aligned'),
    'sigmaa_info': ('aligned'),
    'diff_info': ('aligned'),
    'event_info': ('aligned'),
    'trans_matrix_info': ('aligned'),
    'sdf_info': ('aligned'),
    'single_sdf_file': (''),
    'metadata_info': (''),
    'smiles_info': (''),
    'extra_files': ('extra_files'),
}

# Dictionary containing all references needed to create the zip file
# NB you may need to add a version number to this at some point...
zip_template = {'proteins': {'pdb_info': {},
                             'bound_info': {},
                             'cif_info': {},
                             'mtz_info': {},
                             'diff_info': {},
                             'event_info': {},
                             'sigmaa_info': {},
                             'trans_matrix_info': {}
                             },
                'molecules': {'sdf_files': {},
                              'sdf_info': False,
                              'single_sdf_file': False,
                              'smiles_info': {}
                              },
                'metadata_info': None, }


def _clean_filename(filepath):
    file_split = os.path.splitext(os.path.basename(filepath))
    if fnmatch.fnmatch(
            file_split[0],
            '*_[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'+
            '[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'):
        cleaned_filename = file_split[0][:-8] + file_split[1]
    else:
        cleaned_filename = os.path.basename(filepath)
    return cleaned_filename


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
        cleaned_filename = _clean_filename(filepath)
        ziparchive.write(fullpath,
                         os.path.join(_ZIP_FILEPATHS[param],
                                      cleaned_filename))
        return False

    return True


def _add_file_to_zip_aligned(ziparchive, code, filepath):
    """Add the requested file to the zip archive.

    Args:
        ziparchive: Handle of zip archive
        code: original protein code stripped of any alternate name.
        filepath: filepath from record

    Returns:
        [boolean]: [True of record added to error file]
    """
    media_root = settings.MEDIA_ROOT
    if not filepath:
        return False

    fullpath = os.path.join(media_root, filepath)

    if os.path.isfile(fullpath):
        cleaned_filename = _clean_filename(filepath)
        ziparchive.write(fullpath,
                         os.path.join('aligned', code, cleaned_filename))
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
        with open(combined_sdf_file, 'a') as f_out:
            with open(fullpath, 'r') as f_in:
                f_out.write(f_in.read())
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

        for prot, file in files.items():
            if _add_file_to_zip_aligned(ziparchive, prot.split(":")[0], file):
                error_file.write(
                    '{},{},{}\n'.format(prot, param, file))
                prot_errors += 1
    return prot_errors


def _molecule_files_zip(zip_contents, ziparchive, combined_sdf_file, error_file):
    """Write molecule related data to the ZIP file
       Returns molecule errors
    """

    mol_errors = 0
    for file, prot in zip_contents['molecules']['sdf_files'].items():

        if zip_contents['molecules']['sdf_info'] is True:
            # Add sdf file on the Molecule record to the archive folder.
            if _add_file_to_zip_aligned(ziparchive, prot.split(":")[0], file):
                error_file.write(
                    '{},{},{}\n'.format(prot, 'sdf_info', file))
                mol_errors += 1

        # Append sdf file on the Molecule record to the combined_sdf_file.
        if zip_contents['molecules']['single_sdf_file'] is True:
            if _add_file_to_sdf(combined_sdf_file, file):
                error_file.write(
                    '{},{},{}\n'.format(prot, 'single_sdf_file',
                                        file))
                mol_errors += 1

    return mol_errors


def _smiles_files_zip(zip_contents, ziparchive, download_path):
    """Create and write the smiles file to the ZIP file
    """
    smiles_filename = os.path.join(download_path, 'smiles.smi')
    with open(smiles_filename, 'w') as smilesfile:
        for smi in zip_contents['molecules']['smiles_info']:
            smilesfile.write(smi + ',')
    ziparchive.write(
        smiles_filename,
        os.path.join(_ZIP_FILEPATHS['smiles_info'],
                     os.path.basename(smiles_filename)))
    os.remove(smiles_filename)


def _extra_files_zip(ziparchive, target):
    """If an extra info folder exists at the target root level, then
    copy the contents to the output file as is.
    Note that this will always be the latest information - even for
    preserved searches.
    """

    extra_files = os.path.join(settings.MEDIA_ROOT, 'targets', target.title,
                              'extra_files')
    if os.path.isdir(extra_files):
        for dirpath, dirs, files in os.walk(extra_files):
            for file in files:
                filepath = os.path.join(dirpath, file)
                ziparchive.write(
                    filepath,
                    os.path.join(_ZIP_FILEPATHS['extra_files'], file))


def _create_structures_zip(target,
                           zip_contents,
                           file_url):
    """Write a ZIP file containing data from an input dictionary

    Args:
        target
        zip_contents
        file_url

    """

    logger.info('+ _create_structures_zip')
    download_path = os.path.dirname(file_url)
    os.makedirs(download_path, exist_ok=True)

    error_filename = os.path.join(download_path, "errors.csv")
    error_file = open(error_filename, "w")
    error_file.write("Param, Code, Invalid file reference\n")
    errors = 0

    # If a single sdf file is also wanted then create file to
    # add sdf files to a file called {targer}_combined.sdf.
    combined_sdf_file = None
    if zip_contents['molecules']['single_sdf_file'] is True:
        combined_sdf_file = \
            os.path.join(download_path,
                         '{}_combined.sdf'.format(target.title))

    with zipfile.ZipFile(file_url, 'w') as ziparchive:
        # Read through zip_contents to compile the file
        errors += _protein_files_zip(zip_contents, ziparchive, error_file)

        if zip_contents['molecules']['sdf_files']:
            errors += _molecule_files_zip(zip_contents, ziparchive,
                                          combined_sdf_file, error_file)

        # Add combined_sdf_file to the archive.
        if zip_contents['molecules']['single_sdf_file'] is True \
                and os.path.isfile(combined_sdf_file) :
            ziparchive.write(
                combined_sdf_file,
                os.path.join(_ZIP_FILEPATHS['single_sdf_file'],
                             os.path.basename(combined_sdf_file)))
            os.remove(combined_sdf_file)

        # If smiles info is required, then write one column for each molecule
        # to a smiles.smi file and then add to the archive.
        if zip_contents['molecules']['smiles_info']:
            _smiles_files_zip(zip_contents, ziparchive, download_path)

        # Add the metadata file from the target
        if zip_contents['metadata_info']:
            if _add_file_to_zip(ziparchive, 'metadata_info',
                                zip_contents['metadata_info']):
                error_file.write(
                    '{},{},{}\n'.format(target, 'metadata_info',
                                        zip_contents['metadata_info']))
                errors += 1

        _extra_files_zip(ziparchive, target)

        error_file.close()
        if errors > 0:
            ziparchive.write(error_filename, 'errors.csv')
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
    return [p for p in proteins
            if not fnmatch.fnmatch(p['code'], 'references_*')]


def _create_structures_dict(target, proteins, protein_params, other_params):
    """Write a ZIP file containing data from an input dictionary

    Args:
        target
        proteins
        protein_params
        other_params

    Returns:
        [dict]: [dictionary containing the file contents]
    """
    zip_contents =  copy.deepcopy(zip_template)

    logger.info('+ _create_structures_dict')

    molecules = Molecule.objects.filter(prot_id__in=[protein['id'] for protein
                                                     in proteins]).values()

    # Read through zip_params to compile the parameters
    for protein in proteins:
        for param in protein_params:
            if protein_params[param] is True:
                zip_contents['proteins'][param].update({protein['code']: protein[param]})

    if other_params['single_sdf_file'] is True:
        zip_contents['molecules']['single_sdf_file'] = True

    if other_params['sdf_info'] is True:
        zip_contents['molecules']['sdf_info'] = True

    # sdf information is held as a file on the Molecule record.
    if other_params['sdf_info'] is True or \
            other_params['single_sdf_file'] is True:
        for molecule in molecules:
            protein = Protein.objects.get(id=molecule['prot_id_id'])
            zip_contents['molecules']['sdf_files'].update({molecule['sdf_file']: protein.code})

    # The smiles at molecule level may not be unique.
    if other_params['smiles_info'] is True:
        for molecule in molecules:
            zip_contents['molecules']['smiles_info'].update({molecule['smiles']: None})

    # Add the metadata file from the target
    if other_params['metadata_info'] is True:
        zip_contents['metadata_info'] = target.metadata.name

    return zip_contents


def get_keep_until():
    """Return the keep_until time"""
    return datetime.datetime.now(datetime.timezone.utc) + \
           datetime.timedelta(hours=1)


def get_download_params(request):
    """Check whether structures have been previously downloaded

    Args:
        request

    Returns:
        protein_params, other_params
    """

    protein_param_flags = ['pdb_info', 'bound_info',
                           'cif_info', 'mtz_info',
                           'diff_info', 'event_info',
                           'sigmaa_info', 'trans_matrix_info']
    other_param_flags = ['sdf_info', 'single_sdf_file',
                         'metadata_info', 'smiles_info']

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

    return protein_params, other_params


def check_download_links(request,
                         target,
                         proteins):
    """Check/create the download zip file for dynamic links

    Args:
        request
        target
        proteins

    Returns:
        [file]: [URL to the file in the media directory]
    """


    logger.info('+ _check_download_links')

    protein_params, other_params = get_download_params(request)

    # Remove 'references_' from protein list if present.
    proteins = _protein_garbage_filter(proteins)

    # Save the list of protein codes - this is the ispybsafe set for
    # this user.
    proteins_list = [p['code'] for p in proteins]

    # For dynamic files, if the target, proteins, protein_params and other_params
    # are in an existing record then this is a repeat search .
    # If the zip is there, then we can just return the file_url.
    # If the record is there but the zip file isn't, then regenerate the zip file.
    # from the search, to contain the latest information.
    # If no record is not there at all, then this is a new link.

    if 'static_link' in request.data:
        if request.data['static_link'] is True or \
                request.data['static_link'] == 'true':
            static_link = True
        else:
            static_link = False
    else:
        static_link = False

    existing_link = DownloadLinks.objects.filter(target_id=target.id)\
        .filter(proteins=proteins_list) \
        .filter(protein_params=protein_params)\
        .filter(other_params=other_params) \
        .filter(static_link=False)

    if existing_link:
        if (existing_link[0].zip_file
                and os.path.isfile(existing_link[0].file_url)
                and not static_link):
            return existing_link[0].file_url, True
        elif (os.path.isfile(existing_link[0].file_url) \
                and not static_link):
            # Repeat call, file is currently being created by another process.
            return existing_link[0].file_url, False
        else:
            # Recreate file.
            zip_contents = _create_structures_dict(target,
                                                   proteins,
                                                   protein_params,
                                                   other_params)
            _create_structures_zip(target,
                                   zip_contents,
                                   existing_link[0].file_url)
            existing_link[0].keep_zip_until = get_keep_until()
            existing_link[0].zip_file = True
            if static_link:
                # Changing from a dynamic to a static link
                existing_link[0].static_link = static_link
                existing_link[0].zip_contents = zip_contents
            existing_link[0].save()
            return existing_link[0].file_url, True

    # Create a new link for this user
    # /code/media/downloads/<download_uuid>
    filename = target.title +'.zip'
    media_root = settings.MEDIA_ROOT
    file_url = os.path.join(media_root, 'downloads', str(uuid.uuid4()), filename)

    zip_contents = _create_structures_dict(target,
                                           proteins,
                                           protein_params,
                                           other_params)

    _create_structures_zip(target,
                           zip_contents,
                           file_url)

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
    download_link.save()

    return file_url, True

def recreate_static_file (existing_link):
    """Recreate static file for existing link
    """
    target = existing_link.target

    _create_structures_zip(target,
                           existing_link.zip_contents,
                           existing_link.file_url)
    existing_link.keep_zip_until = get_keep_until()
    existing_link.zip_file = True
    existing_link.save()


def maintain_download_links():
    """Maintain zip files

       Physical zip files are removed after 1 hour (or when the next POST call
       is made) for security reasons and to conserve memory space.

    """

    old_links = DownloadLinks.objects.\
        filter(keep_zip_until__lt=datetime.datetime.now(datetime.timezone.utc)).\
        filter(zip_file=True)
    for old_link in old_links:
        old_link.zip_file = False
        old_link.save()
        shutil.rmtree(os.path.dirname(old_link.file_url), ignore_errors=True)
