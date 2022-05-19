"""Functions (for the celery tasks) that transfer (download) files from
Fragalysis to Squonk.
"""
import os
import shutil

from django.conf import settings
from rest_framework import status
from dm_api.dm_api import DmApi

from rdkit import Chem

from celery.utils.log import get_task_logger
from viewer.utils import (
    clean_filename,
    add_prop_to_sdf,
    add_prop_to_mol,
    create_media_sub_directory,
    delete_media_sub_directory
)
from viewer.models import (
    Molecule,
    Protein,
    ComputedMolecule,
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

SQUONK_PROT_MAPPING = {
    'apo_desolve': {'func': 'prot_file', 'field': 'apo_desolve_info'},
    'mol': {'func': 'mol_field', 'field': 'sdf_info'},
    'sdf': {'func': 'sdf_file', 'field': 'sdf_file'}
}

SQUONK_COMP_MAPPING = {
    'mol': {'func': 'comp_mol_field', 'field': 'sdf_info'},
}

REF_PROP = 'ref_mols'

def ref_mol_from_protein_code(protein_code, target_title):
    """Returns the xchem ref_mol given a protein code and target name
    If it can't find a target prefix it just returns the code as it is.

    Args:
        protein_code (e.g. Mpro-x11513_0A:ALP-POS-f13221e1-4 )
    Returns:
        xchem ref mol (e.g. x11513_0A)
    """
    target_prefix = target_title + '-'
    if target_prefix in protein_code:
        return protein_code.strip().split(":")[0].replace(target_prefix, '')
    else:
        return protein_code


def mol_field(field, trans_dir, protein_code, target):
    """Extract the contents of the mol_field to a file in the transfer directory

    Args:
        field : source field in Molecule table
        trans_dir : temporary directory used for transfer
        protein
        target
    Returns:
        filepath
    """
    protein = Protein.objects.get(code=protein_code)
    mol = Molecule.objects.get(prot_id=protein.id)

    mol_block = getattr(mol, field)
    filepath = os.path.join(trans_dir, protein.code.strip().split(":")[0] + '.mol')
    code = ref_mol_from_protein_code(protein_code, target.title)
    add_prop_to_mol(mol_block, filepath, code)

    return filepath


def sdf_file(field, trans_dir, protein_code, target):
    """Copy the requested file from the molecule table to the transfer directory, standardising
    the filename.

    Args:
        field : source field in Molecule table
        trans_dir : temporary directory used for transfer
        protein
        target
    Returns:
        filepath
    """

    protein = Protein.objects.get(code=protein_code)
    mol = Molecule.objects.get(prot_id=protein.id)
    file = getattr(mol, field)
    if not file:
        return None

    inpath = os.path.join(settings.MEDIA_ROOT, file.name)
    filepath = os.path.join(trans_dir, clean_filename(inpath))
    code = ref_mol_from_protein_code(protein_code, target.title)
    add_prop_to_sdf(inpath, filepath, {REF_PROP: code})
    return filepath


def prot_file(field, trans_dir, protein_code, target):
    """Copy the requested file from the protein table to the transfer directory, standardising the
    name.

    Args:
        field : source field in Protein table
        trans_dir : temporary directory used for transfer
        protein
        target
    Returns:
        filepath
    """

    protein = Protein.objects.get(code=protein_code)
    file = getattr(protein, field)
    if not file:
        return None

    inpath = os.path.join(settings.MEDIA_ROOT, file.name)
    filepath = os.path.join(trans_dir, clean_filename(inpath))
    shutil.copyfile(inpath, filepath)
    return filepath


def comp_mol_field(field, trans_dir, name, target):
    """Extract the contents of the mol_field to a file in the transfer directory

    Args:
        field : source field in Molecule table
        trans_dir : temporary directory used for transfer
        protein
        target
    Returns:
        filepath
    """
    comp = ComputedMolecule.objects.get(name=name)
    mol_block = getattr(comp, field)
    filepath = os.path.join(trans_dir, name + '.mol')
    # In the case of a computed molecule the whole name is used.
    add_prop_to_mol(mol_block, filepath, name)
    return filepath


def process_file_transfer(auth_token,
                          job_transfer_id):
    """Check/create the file transfer for list of proteins

    Args:
        request,
        job_transfer_id

    """

    logger.info('+ process_file_transfer')
    logger.debug(auth_token)
    job_transfer = JobFileTransfer.objects.get(id=job_transfer_id)

    trans_sub_dir = os.path.join('squonk_transfer', str(job_transfer.id))
    trans_dir = create_media_sub_directory(trans_sub_dir)

    # location in squonk project where files will reside e.g. "fragalysis-files"
    target = job_transfer.target
    squonk_directory = settings.SQUONK_MEDIA_DIRECTORY + '/' + target.title

    # This to pick up NULL values from the changeover to using compounds.
    if not job_transfer.compounds:
        job_transfer.compounds = []

    num_to_transfer = len(job_transfer.proteins) + len(job_transfer.compounds)
    idx = 0

    # Proteins
    for protein_code in job_transfer.proteins:
        # For each protein transfer the list of files to squonk
        # Then update the progress in the job_transfer record
        logger.info('+ Processing files for: ' + protein_code)
        file_list = []
        for file_type in SQUONK_PROT_MAPPING.values():
            filepath = globals()[file_type['func']](file_type['field'], trans_dir,
                                                          protein_code, target)
            if filepath:
                file_list.append(filepath)

        logger.info(squonk_directory)
        logger.info(file_list)

        result = DmApi.put_unmanaged_project_files(access_token=auth_token,
                                                   project_id=job_transfer.squonk_project,
                                                   project_files=file_list,
                                                   project_path=squonk_directory,
                                                   force=True)
        logger.info(result)

        if result.success:
            idx += 1
            job_transfer.transfer_progress = (idx*100/num_to_transfer)
            job_transfer.save()
        else:
            logger.info('File Transfer Failed')
            raise RuntimeError('File Transfer Failed')


    # Computed Molecules
    for name in job_transfer.compounds:
        # For each protein transfer the list of files to squonk
        # Then update the progress in the job_transfer record
        logger.info('+ Processing files for: ' + name)
        file_list = []
        for file_type in SQUONK_COMP_MAPPING.values():
            filepath = globals()[file_type['func']](file_type['field'], trans_dir,
                                                          name, target)
            if filepath:
                file_list.append(filepath)

        logger.info(squonk_directory)
        logger.info(file_list)

        result = DmApi.put_unmanaged_project_files(access_token=auth_token,
                                                   project_id=job_transfer.squonk_project,
                                                   project_files=file_list,
                                                   project_path=squonk_directory,
                                                   force=True)
        logger.info(result)

        if result.success:
            idx += 1
            job_transfer.transfer_progress = (idx*100/num_to_transfer)
            job_transfer.save()
        else:
            logger.error('File Transfer Failed')
            raise RuntimeError('File Transfer Failed')

    # Tidy up the transfer directory.
    delete_media_sub_directory(trans_sub_dir)


def check_file_transfer(request):
    """Check the request and return a list of protein codes and/or computed molecule names

    Args:
        request
    Returns
        error dict
        dict of validated protein codes
        dict of validated computed molecule names
    """
    error  = {}
    proteins = []
    compounds = []

    if request.data['proteins']:
        # Get first part of protein code
        proteins_list = [p.strip().split(":")[0]
                         for p in request.data['proteins'].split(',')]
        proteins = []
        # Filter by protein codes
        for code_first_part in proteins_list:
            prot = Protein.objects.filter(code__contains=code_first_part).values()
            if prot.exists():
                proteins.append(prot.first())
            else:
                error['message'] = 'Please enter valid protein code for' \
                                   + ': {} '.format(code_first_part)
                error['status'] = status.HTTP_404_NOT_FOUND
                return error, proteins, compounds

        if len(proteins) == 0:
            error['message'] = 'API expects a list of comma-separated protein codes'
            error['status'] = status.HTTP_404_NOT_FOUND
            return error, proteins, compounds

    if request.data['compounds']:
        # Get first part of protein code
        compounds_list = [c.strip() for c in request.data['compounds'].split(',')]
        compounds = []
        # Filter by protein codes
        for compound in compounds_list:
            comp = ComputedMolecule.objects.filter(name=compound).values()
            if comp.exists():
                compounds.append(comp.first())
            else:
                error['message'] = 'Please enter valid compound name for' \
                                   + ': {} '.format(compound)
                error['status'] = status.HTTP_404_NOT_FOUND
                return error, proteins, compounds

        if len(compounds) == 0:
            error['message'] = 'API expects a list of comma-separated compound names'
            error['status'] = status.HTTP_404_NOT_FOUND
            return error, proteins, compounds

    if proteins or compounds:
        return error, proteins, compounds
    else:
        error['message'] = 'A valid protein codes and/or a list of valid compound names must ' \
                           'be entered'
        error['status'] = status.HTTP_404_NOT_FOUND
        return error, proteins, compounds
