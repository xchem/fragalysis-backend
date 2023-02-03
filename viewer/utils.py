"""
utils.py

Collection of technical methods tidied up in one location.
"""
import fnmatch
import os
import shutil

from django.conf import settings

from rdkit import Chem

# Set .sdf file format version
# Used at the start of every SDF file.
SDF_VERSION = 'ver_1.2'


def create_squonk_job_request_url(instance_id):
    """Creates the Squonk Instance API url from an instance ID (UUID).
    """
    return settings.SQUONK2_INSTANCE_API + str(instance_id)


def create_media_sub_directory(sub_directory):
    """Creates a directory (or directories) in the MEDIA directory,
    returning the full path.
    """
    directory = os.path.join(settings.MEDIA_ROOT, sub_directory)
    os.makedirs(directory, exist_ok=True)

    return directory


def delete_media_sub_directory(sub_directory):
    """Removes a media sub-directory."""
    assert sub_directory
    assert len(sub_directory)
    directory = os.path.normpath(os.path.join(settings.MEDIA_ROOT, sub_directory))
    if not os.path.isdir(directory):
        # No such directory!
        return
    if not directory.startswith(settings.MEDIA_ROOT) or \
            os.path.samefile(directory, settings.MEDIA_ROOT):
        # Danger!
        return

    shutil.rmtree(directory)


def add_prop_to_sdf(sdf_file_in, sdf_file_out, properties):
    """Returns an SDF file with the requested property
    (a dictionary of keys and values) added.
    Note that this assumes the file starts with a blank molecule.

    SDF parameters are of format:
        >  <TransFSScore>  (1)
        0.115601
    In this case we'd provide the following properties dictionary...
        {"TransFSScore": "0.115601"}

    Args:
        sdf_file_in : filepath of source file
        sdf_file_out : filepath to write new file to
        properties : dictionary or properties
    """
    _REC_SEPARATOR = '$$$$\n'

    found_separator = False
    with open(sdf_file_out, 'a', encoding='utf-8') as sdf_out:
        with open(sdf_file_in, 'r', encoding='utf-8') as sdf_in:
            while True:
                line = sdf_in.readline()
                if line:
                    if not found_separator and line == _REC_SEPARATOR:
                        # Found first separator
                        # dump the parameters now
                        found_separator = True
                        for name, value in properties.items():
                            sdf_out.write(f'>  <{name}>  (1)\n')
                            sdf_out.write(f'{value}\n')
                            sdf_out.write('\n')
                        # Follow with separator
                        sdf_out.write(line)
                    else:
                        sdf_out.write(line)
                else:
                    break


def add_prop_to_mol(mol_field, mol_file_out, value):
    """Returns a mol_file with the requested property added
    Note that the only property that seems to work with a .mol file if you write it is: "_Name".

    Args:
        mol_field
        mol_file_out : filepath to write new file to
        property
        value
    """
    rd_mol = Chem.MolFromMolBlock(mol_field)
    # Set the new property in the property mol object
    rd_mol.SetProp("_Name", value)
    Chem.MolToMolFile(rd_mol, mol_file_out)


def clean_filename(filepath):
    """Return the "clean" version of a Django filename without the '_abcdefg_' that is created
    when a file is overwritten.

    Given a path and a file, e.g. './media/sdfs/Mpro-x3351_0A_rtEVbqf.sdf' this
    function will return 'Mpro-x3351_0A.sdf'.

    Args:
        filepath

    Returns:
        cleaned filename
    """
    file_split = os.path.splitext(os.path.basename(filepath))
    if fnmatch.fnmatch(
            file_split[0],
            '*_[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'+
            '[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'):
        cleaned_filename = file_split[0][:-8] + file_split[1]
    else:
        cleaned_filename = os.path.basename(filepath)
    return cleaned_filename


def get_https_host(request):
    """Common enabler code for returning https urls

    This is to map links to HTTPS to avoid Mixed Content warnings from Chrome browsers
    SECURE_PROXY_SSL_HEADER is referenced because it is used in redirecting URLs - if
    it is changed it may affect this code.
    Using relative links will probably also work, but This workaround allows both the
    'download structures' button and the DRF API call to work.
    Note that this link will not work on local
    """
    return settings.SECURE_PROXY_SSL_HEADER[1] + '://' + request.get_host()
