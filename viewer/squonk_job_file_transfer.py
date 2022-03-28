"""
squonk_job_file_transfer Functions for the celery tasks for transferring
files from Fragalysis to Squonk.

"""
import os
import shutil

from django.conf import settings

from dm_api.dm_api import DmApi

from celery.utils.log import get_task_logger
from viewer.utils import clean_filename

from viewer.models import (
    Molecule,
    Protein,
    JobFileTransfer
)
logger = get_task_logger(__name__)

# Squonk Specification File Paths
# squonk-file-paths:
#   lhs-fragments:
#     molecule: {pdb_code}.mol
#     template: /fragalysis-files/{target_name}/{pdb_code}.mol
#   lhs-protein-apo-desolv:
#     molecule: {pdb_code}_apo_desolv.pdb
#     template: /fragalysis-files/{target_name}/{pdb_code}_apo_desolv.pdb

# Fragalysis Upload Mapping (as well as I can decipher it:
#
# Type         | From Aligned File     |  E.g.                          | Stored In
# .apo_desolve | <pcode>_apo-desolv.pdb| CD44MMA-x0017_0A_apo-desolv.pdb| Protein.apo_desolve_info
# .mol         | <pcode>.mol           | CD44MMA-x0017_0A.mol           | Molecule.sdf_info
# .sdf         | <pcode>.sdf           | CD44MMA-x0017_0A.sdf           | Molecule.sdf_file
#

SQUONK_MAPPING = {
    'apo_desolve': { 'func': 'prot_file', 'field': 'apo_desolve_info'},
    'mol': { 'func': 'mol_field', 'field': 'sdf_info'},
    'sdf': { 'func': 'mol_file', 'field': 'sdf_file'}
}

def mol_field(field, trans_dir, protein_code):
    """Extract the contents of the mol_field to a file in the transfer directory

    Args:
        field : source field in Molecule table
        trans_dir : temporary directory used for transfer
        protein
    Returns:
        filepath
    """
    protein = Protein.objects.get(code=protein_code)
    mol = Molecule.objects.get(prot_id=protein.id)
    file = getattr(mol, field)
    filepath = os.path.join(trans_dir, protein.code.strip().split(":")[0] + '.mol')
    with open(filepath, 'a') as f_out:
        f_out.write(file)
    return filepath


def mol_file(field, trans_dir, protein_code):
    """Copy the requested file from the molecule table to the transfer directory, standardising
    the filename.

    Args:
        field : source field in Molecule table
        trans_dir : temporary directory used for transfer
        protein
    Returns:
        filepath
    """

    protein = Protein.objects.get(code=protein_code)
    mol = Molecule.objects.get(prot_id=protein.id)
    file = getattr(mol, field)
    inpath = os.path.join(settings.MEDIA_ROOT, file.name)
    filepath = os.path.join(trans_dir, clean_filename(inpath))
    shutil.copyfile(inpath, filepath)
    return filepath


def prot_file(field, trans_dir, protein_code):
    """Copy the requested file from the protein table to the transfer directory, standardising the
    name.

    Args:
        field : source field in Protein table
        trans_dir : temporary directory used for transfer
        protein
    Returns:
        filepath
    """

    protein = Protein.objects.get(code=protein_code)
    file = getattr(protein, field)
    inpath = os.path.join(settings.MEDIA_ROOT, file.name)
    filepath = os.path.join(trans_dir, clean_filename(inpath))
    shutil.copyfile(inpath, filepath)
    return filepath


def process_file_transfer(auth_token,
                          job_transfer_id):
    """Check/create the file transfer for list of proteins

    Args:
        request,
        job_transfer_id

    """

    logger.info('+ process_file_transfer')
    logger.info (auth_token)
    job_transfer = JobFileTransfer.objects.get(id=job_transfer_id)
    #job_transfer.transfer_spec = list(SQUONK_MAPPING.keys())
    #job_transfer.save()

    trans_dir = os.path.join(settings.MEDIA_ROOT, 'squonk_transfer', str(job_transfer.id))
    target = job_transfer.target
    # Format e.g.: /fragalysis-files/Mpro
    squonk_directory = settings.SQUONK_MEDIA_DIRECTORY + target.title + '/'
    os.makedirs(trans_dir, exist_ok=True)

    len_proteins = len(job_transfer.proteins)
    idx = 0

    for protein_code in job_transfer.proteins:
        # For each protein transfer the list of files to squonk
        # Then update the progress in the job_transfer record
        logger.info('+ Processing files for: ' + protein_code)
        file_list = []
        for file_type in SQUONK_MAPPING.values():
            file_list.append(globals()[file_type['func']](file_type['field'], trans_dir,
                                                          protein_code))

        logger.info(squonk_directory)
        logger.info(file_list)

        result = DmApi.upload_unmanaged_project_files(access_token=auth_token,
                                                      project_id=job_transfer.squonk_project,
                                                      project_files=file_list)
#                                                     project_path=squonk_directory,
#                                                     force=True)
        logger.info(result)

        if result.success:
            idx += 1
            job_transfer.transfer_progress = (idx*100/len_proteins)
            job_transfer.save()
        else:
            logger.info('File Transfer Failed')
            raise RuntimeError('File Transfer Failed')

    # Tidy up the transfer directory.
    shutil.rmtree(trans_dir)
