import datetime
import os
import psutil
import shutil

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")

import django
django.setup()
from django.conf import settings
from celery import shared_task
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors
from viewer.models import Target, Compound, DesignSet
import zipfile

from .sdf_check import (
    add_warning,
    check_SMILES,
    check_ver_name,
    check_blank_mol_props,
    check_blank_prop,
    check_mol_props,
    check_name_characters,
    check_refmol,
    check_field_populated,
    check_compound_set
)
from .compound_set_upload import get_additional_mols
from .target_set_upload import process_target, validate_target
from .cset_upload import blank_mol_vals, MolOps, PdbOps
from .squonk_job_file_transfer import process_file_transfer, SQUONK_MAPPING
from .models import ComputedSet, JobFileTransfer, Molecule, User

from celery.utils.log import get_task_logger
logger = get_task_logger(__name__)


def check_services():
    """ Method to ensure redis and celery services are running to allow handling of tasks - attempts to start either
    service if it is not running, and returns True or False to indicate whether both services are running

    Parameters
    ----------
    None

    Returns
    -------
    status: bool
        True if both services are running, False if one or both are not

    """
    services = [p.name() for p in psutil.process_iter()]
    if 'redis-server' not in services:
        os.system('redis-server &')
    if 'celery' not in services:
        os.system('celery -A fragalysis worker -l info &')
    services = [p.name() for p in psutil.process_iter()]
    if 'redis-server' not in services or 'celery' not in services:
        return False
    return True


# Uploading Compound Sets ###

@shared_task
def process_compound_set(validate_output):
    """Celery task to process a computed set. Used in a chain invoked from views.py's
    UploadCSet, taking the output of the `validate_compound_set()`, and uploads
    molecules to a new computed set if the uploaded files are valid.

    (Code is legacy, and has been taken on by multiple users - comments often
    retro-fitted and no apologies for the design)

    Parameters
    ----------
    validate_output: tuple
        see return value from 'validate_compound_set()'

    Returns
    -------
    If successful

        process_stage: str
            'process'
        process_type: str
            'cset'
        compound_set.name: str
            name of the computed set

    Otherwise (i.e. not validated)...

        Returns the response from 'validate_compound_set()'
    """

    process_stage, process_type, validate_dict, validated, params = validate_output

    logger.info('process_compound_set() ENTER')
    logger.info('process_compound_set() validated=%s', validated)
    logger.info('process_compound_set() params=%s', params)

    if not validated:
        logger.warning('process_compound_set() EXIT (not validated)')
        return process_stage, 'cset', validate_dict, validated


    submitter_name, submitter_method, blank_version = blank_mol_vals(params['sdf'])
    zfile, zfile_hashvals = PdbOps().run(params)

    # This used to be in 'validate' (incorrectly), and is now here.
    # Create a unique name and remove any that existing set with that name.
    # Basically uploading erases any prior set.
    unique_name = "".join(submitter_name.split()) + '-' + "".join(submitter_method.split())
    csets = ComputedSet.objects.filter(unique_name=unique_name)
    for cset in csets:
        cset.delete()

    save_mols = MolOps(user_id=params['user_id'],
                       sdf_filename=params['sdf'],
                       submitter_name=submitter_name,
                       submitter_method=submitter_method,
                       target=params['target'],
                       version=blank_version,
                       zfile=zfile,
                       zfile_hashvals=zfile_hashvals)
    compound_set = save_mols.task()

    logger.warning('process_compound_set() EXIT (%s)', compound_set.name)
    return 'process', 'cset', compound_set.name

# End Uploading Compound Sets ###

# Validating Compound Sets ###

# Set .sdf format version here
version = 'ver_1.2'


@shared_task
def validate_compound_set(user_id, sdf_file, target=None, zfile=None, update=None):
    """ Celery task to process validate the uploaded files for a computed set upload.
    SDF file is mandatory, zip file is optional. Used in a chain invoked from views.py's
    UploadCSet and the first stage in a chain that includes 'process_compound_set()'

    (Code is legacy, and has been taken on by multiple users - comments often
    retro-fitted and no apologies for the design)

    Parameters
    ----------
    user_id: The record ID of the user initiating the task.
    sdf_file: str
        filepath of the uploaded sdf file, which is saved to temporary storage
        by `viewer.views.UploadCSet`
    target: str
        name of the target (`viewer.models.Target.title`) to add add the computed set to
    zfile: dict
        dictionary where key is the name of the file minus extension and path,
        and value is the filename, which is saved
        to temporary storage by `viewer.views.UploadCSet`
    update: str
        The name of a computed set to update, or 'None'

    Returns
    -------
    tuple containing the following:
        - processing stage (str): 'validate'
        - processing type (str): 'cset'
        - validate dict (dict): dict containing any errors found during the calidation step
        - validated (bool): True if the file(s) were validated, False if not
        - params (dict):
            - user_id (int): User record ID of user initiating the task
            - sdf (str): name of the uploaded sdf file
            - target (str): name of the target that the computed set is associated with
            - zfile: dictionary where key is the name of the file minus extension
                     and path, and value is the filename, which is saved to
                     temporary storage by `viewer.views.UploadCSet`
    If there are problems with the 'blank (initial) molecule' the following is also
    returned...
            - submitter_name (str): name of the author of the computed set
            - submitter_method (str): name of the method used to generate the computed set
    """

    logger.info('validate_compound_set() ENTER')
    logger.info('validate_compound_set() user_id=%s', user_id)
    logger.info('validate_compound_set() sdf_file=%s', sdf_file)
    logger.info('validate_compound_set() target=%s', target)
    logger.info('validate_compound_set() update=%s', update)

    validated = True
    validate_dict = {'molecule_name': [],
                     'field': [],
                     'warning_string': []}

    suppl = Chem.SDMolSupplier(sdf_file)
    # print('%d mols detected (including blank mol)' % (len(suppl),))
    blank_mol = suppl[0]

    # Get submitter name/info for passing into upload to get unique name
    submitter_name = blank_mol.GetProp('submitter_name')
    submitter_method = blank_mol.GetProp('method')

    if blank_mol is None:
        validate_dict = add_warning(molecule_name='Blank Mol',
                                    field='N/A',
                                    warning_string='your blank molecule could not be'
                                                   ' read by rdkit. The molecule must'
                                                   ' have at least one atom!'
                                                   ' No other checks were done',
                                    validate_dict=validate_dict)
        validated = False
        logger.info('validate_compound_set() EXIT validated=%s', validated)
        return (validate_dict, validated, sdf_file, target, zfile,
                submitter_name, submitter_method)

    if not update or update == 'None':
        validate_dict = check_compound_set(blank_mol, validate_dict)
    else:
        validate_dict = check_compound_set(blank_mol, validate_dict, update=update)

    other_mols = []
    for i in range(1, len(suppl)):
        other_mols.append(suppl[i])

    # all mol checks
    # Check if all mols can be read by rdkit
    # - all mols have the same properties
    all_props = []
    # Use index and check_mol to see if any sdf entries are None type mol objects
    index = 1
    for mol in suppl:
        if not mol:
            add_warning(molecule_name='Unknown',
                        field='N/A',
                        warning_string=f'SDF entry number: {index} cannot be converted'
                                        ' into an rdkit mol object',
                        validate_dict=validate_dict)
        if mol:
            all_props.extend([key for key in list(mol.GetPropsAsDict().keys())])
        index += 1
    unique_props = list(set(all_props))

    for mol in suppl:
        if mol:
            props = [key for key in list(mol.GetPropsAsDict().keys())]
            diff_list = np.setdiff1d(props, unique_props)
            for diff in diff_list:
                add_warning(molecule_name=mol.GetProp('_Name'),
                            field='property (missing)',
                            warning_string=f'{diff} property is missing from this molecule',
                            validate_dict=validate_dict)

    # Check version in blank mol
    validate_dict = check_ver_name(blank_mol, version, validate_dict)

    # Check compulsory fields in blank mol props
    validate_dict = check_blank_mol_props(blank_mol, validate_dict)

    # Check properties have been described and validate url
    # Anything wrong with validation goes into this and we leave
    # with a validation error.
    validate_dict = check_blank_prop(blank_mol, validate_dict)

    # main mols checks
    # - missing compulsary fields
    # - check name characters
    # - check pdb assignment and if pdb filepath exists
    # - check compulsory field populated
    # - check SMILES can be opened by rdkit
    # (check api for pdb if fragalysis)
    for m in other_mols:
        if m:
            validate_dict = check_mol_props(m, validate_dict)
            validate_dict = check_name_characters(m.GetProp('_Name'), validate_dict)
            # validate_dict = check_pdb(m, validate_dict, target, zfile)
            validate_dict = check_refmol(m, validate_dict, target)
            validate_dict = check_field_populated(m, validate_dict)
            validate_dict = check_SMILES(m, validate_dict)

    len_validate_dict = len(validate_dict['molecule_name'])
    if len_validate_dict != 0:
        logger.warning('validate_compound_set() validated=False'
                       ' len_validate_dict=%s', len_validate_dict)
        validated = False

    # params = {
    #     'sdf': 'tests/test_data/test_chodera/test_0_2.sdf',
    #     'target': 'Mpro',
    #     'choice': 1,
    #     'update_set': None,
    #     'pdb_zip': '/code/tests/test_data/test_chodera/references.zip'
    # }

    params = {
        'user_id': user_id,
        'sdf': sdf_file,
        'target': target,
        'pdb_zip': zfile
    }

    logger.info('validate_compound_set() EXIT validated=%s params=%s', validated, params)
    return 'validate', 'cset', validate_dict, validated, params

# End Validating Compound Sets ###

# Design sets ###

def create_mol(inchi, long_inchi=None, name=None):
    # check for an existing compound
    if long_inchi:
        cpd = Compound.objects.filter(long_inchi=long_inchi)
        sanitized_mol = Chem.MolFromInchi(long_inchi, sanitize=True)
    else:
        cpd = Compound.objects.filter(inchi=inchi)
        sanitized_mol = Chem.MolFromInchi(inchi, sanitize=True)

    if len(cpd) != 0:
        new_mol = cpd[0]
    else:
        # add molecule and return the object
        new_mol = Compound()

    new_mol.smiles = Chem.MolToSmiles(sanitized_mol)
    new_mol.inchi = inchi
    if long_inchi:
        new_mol.long_inchi = long_inchi
    new_mol.identifier = name

    # descriptors
    new_mol.mol_log_p = Chem.Crippen.MolLogP(sanitized_mol)
    new_mol.mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(sanitized_mol))
    new_mol.heavy_atom_count = Chem.Lipinski.HeavyAtomCount(sanitized_mol)
    new_mol.heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(sanitized_mol))
    new_mol.nhoh_count = Chem.Lipinski.NHOHCount(sanitized_mol)
    new_mol.no_count = Chem.Lipinski.NOCount(sanitized_mol)
    new_mol.num_h_acceptors = Chem.Lipinski.NumHAcceptors(sanitized_mol)
    new_mol.num_h_donors = Chem.Lipinski.NumHDonors(sanitized_mol)
    new_mol.num_het_atoms = Chem.Lipinski.NumHeteroatoms(sanitized_mol)
    new_mol.num_rot_bonds = Chem.Lipinski.NumRotatableBonds(sanitized_mol)
    new_mol.num_val_electrons = Descriptors.NumValenceElectrons(sanitized_mol)
    new_mol.ring_count = Chem.Lipinski.RingCount(sanitized_mol)
    new_mol.tpsa = Chem.rdMolDescriptors.CalcTPSA(sanitized_mol)

    # make sure there is an id so inspirations can be added
    new_mol.save()

    return new_mol

def process_design_compound(compound_row):
    # sanitize, generate mol and inchi
    smiles = compound_row['smiles']
    name = compound_row['identifier']
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    sanitized_mol_smiles = Chem.MolToSmiles(mol, canonical=True)
    sanitized_mol = Chem.MolFromSmiles(sanitized_mol_smiles)
    inchi = Chem.inchi.MolToInchi(sanitized_mol)
    long_inchi = None

    if len(inchi)>255:
        # TODO: get_inchi in model
        inchi = str(inchi)[0:255]
        long_inchi = inchi

    new_mol = create_mol(inchi)

    # deal with inspirations
    inspirations = compound_row['inspirations'].split(',')
    for insp in inspirations:
        # TODO: find matching molecules - change to molecules and search history to find the correct version.
        #  -- search all history and find most recent with matching code? or code most closely matching design date?
        # (Can this be accessed, or does the view need changing to display the correct one? Not implemented yet anyway)
        molecules = Molecule.objects.filter(prot_id__code__contains=insp.split(':')[0].split('_')[0])
        # compounds = [m.cmpd_id for m in molecules]
        for molecule in molecules:
            new_mol.inspirations.add(molecule)

    # save the molecule and return it
    new_mol.save()
    return new_mol


def process_one_set(set_df, name, set_type=None, set_description=None):
    # add new set
    new_set = DesignSet()
    new_set.set_name = name
    new_set.set_type = set_type
    new_set.set_description = set_description
    new_set.save()
    compounds = []
    for _, row in set_df.iterrows():
        compounds.append(process_design_compound(row))

    for compound in compounds:
        new_set.compounds.add(compound)

    new_set.save()

    return compounds


def process_design_sets(df, set_type=None, set_description=None):
    # Unused variables
    del set_type
    del set_description

    set_names = list(set(df['set_name']))
    sets = []
    for name in set_names:
        set_df = df[df['set_name'] == name]
        compounds = process_one_set(set_df, name)
        sets.append(compounds)

    return set_names, sets
# End Design sets ###


# Target Sets ###
@shared_task
def validate_target_set(target_zip, target=None, proposal=None, email=None):
    """ Celery task to process validate the uploaded files/format for a target set upload. Zip file is mandatory

    Parameters
    ----------
    target_zip: str
        filepath of the uploaded target file, which is saved to temporary storage by `viewer.views.UploadTSet`
    target: str
        name of the target (`viewer.models.Target.title`) to add add the target set to

    Returns
    -------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded target file
            - target (str): name of the target
            - submitter_name (str): name of the user who submitted the upload

    """
    logger.info('+ validating target set: ' + target_zip)

    # Get submitter name/info for passing into upload to get unique name
    submitter_name = ''

    tmp_folder = os.path.dirname(target_zip)
    new_data_folder = os.path.join(tmp_folder, 'new_data')

    # This will create the target folder in the tmp/ location.
    with zipfile.ZipFile(target_zip, 'r') as zip_ref:
        zip_ref.extractall(new_data_folder)

    validated, validate_dict = validate_target(new_data_folder, target, proposal)

    os.remove(target_zip)
    # Tidy up data if not validated
    if not validated:
        shutil.rmtree(new_data_folder)

    return ('validate', 'tset', validate_dict, validated, new_data_folder, target, proposal,
            submitter_name, email)


@shared_task
def process_target_set(validate_output):
    """ Celery task to process a target set, that takes the output of the validation task, and uploads molecules to a
    new target set if the uploaded files are valid

    Parameters
    ----------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded target file (zip file)
            - target (str): name of the target
            - submitter_name (str): name of the author of the computed set

    Returns
    -------
    target str
        name of the processed target set

    """
    # Validate output is a tuple - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    process_stage, process_type, validate_dict, validated, new_data_folder, target_name, proposal_ref, submitter_name, \
        contact_email = validate_output

    # If there is a validation error, stop here.
    if not validated:
        return process_stage, 'tset', validate_dict, validated

    if validated:
        logger.info('+ processing target set: ' + target_name + ' target_folder:' + new_data_folder)

        mols_loaded, mols_processed = process_target(new_data_folder, target_name, proposal_ref)

        return 'process', 'tset', target_name, mols_loaded, mols_processed, contact_email
# End Target Sets ###


# File Transfer ###
@shared_task
def process_job_file_transfer(auth_token, id):
    """ Celery task to take a list of proteins and specification and transfer the files to Squonk2

    Parameters
    ----------
    id
        id of job_file_transfer record

    Returns
    -------
    target str
        name of the processed target set

    """

    logger.info('+ Starting File Transfer (%s)', id)
    job_transfer = JobFileTransfer.objects.get(id=id)
    job_transfer.transfer_status = "STARTED"
    job_transfer.transfer_task_id = str(process_job_file_transfer.request.id)
    job_transfer.save()
    try:
        process_file_transfer(auth_token, job_transfer.id)
    except RuntimeError as error:
        logger.info('- File Transfer failed %s', id)
        logger.info(error)
        job_transfer.transfer_status = "FAILURE"
        job_transfer.save()
    else:
        logger.info('+ File Transfer Successful (%s)', id)
        # Update the transfer datetime for comparison with the target upload datetime.
        # This should only be done on a successful upload.
        job_transfer.transfer_datetime = datetime.datetime.now(datetime.timezone.utc)
        job_transfer.transfer_progress = 100.00
        job_transfer.transfer_status = "SUCCESS"
        files_spec = list(SQUONK_MAPPING.keys())
        job_transfer.transfer_spec = files_spec
        job_transfer.save()
        logger.info('- File Transfer Ended Successfully (%s)', id)

    return job_transfer.transfer_status

# End File Transfer ###
