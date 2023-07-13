import logging
import datetime
import os
import psutil
import shutil

from django.db import IntegrityError



os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")

import django
django.setup()

from django.conf import settings
from celery import shared_task
from fragalysis.celery import app as celery_app
from celery.utils.log import get_task_logger
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors
from viewer.models import Compound, DesignSet
import zipfile

from viewer.target_loader import load_target


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
from .target_set_upload import process_target, validate_target
from .cset_upload import blank_mol_vals, MolOps, PdbOps
from .squonk_job_file_transfer import (
    process_file_transfer,
    SQUONK_PROT_MAPPING,
    SQUONK_COMP_MAPPING
)
from .squonk_job_file_upload import (
    get_upload_sub_directory,
    process_compound_set_file
)
from .models import ComputedSet, JobRequest, JobFileTransfer, Molecule
from .utils import SDF_VERSION, delete_media_sub_directory

# If Celery configured to always run 'synchronously' (eager),
# then use a standard logger, otherwise use the Celery logger.
if settings.CELERY_TASK_ALWAYS_EAGER:
    logger = logging.getLogger(__name__)
else:
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
    redis_service = 'redis-server'
    if redis_service not in services:
        logger.info('Service "%s" not present ... starting it', redis_service)
        rc = os.system('redis-server &')
        if rc != 0:
            logger.warning('Got exit-code %s starting "%s"', rc, redis_service)

    celery_service = 'celery'
    if settings.CELERY_TASK_ALWAYS_EAGER:
        logger.warning('Found CELERY_TASK_ALWAYS_EAGER - skipping check of celery service')
        return True
    else:
        if celery_service not in services:
            logger.info('Service "%s" not present ... starting it', celery_service)
            rc = os.system('celery -A fragalysis worker -l info &')
            if rc != 0:
                logger.warning('Got exit-code %s starting "%s"', rc, celery_service)

    # Check again...
    services = [p.name() for p in psutil.process_iter()]
    if redis_service not in services:
        logger.error('Redis is not running as a service')
        return False
    if not settings.CELERY_TASK_ALWAYS_EAGER and celery_service not in services:
        logger.error('Celery is not running as a service')
        return False

    # OK if we get here
    logger.info('Redis and Celery are running as services')
    return True


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
    If successful, a dictionary that contains the fields: -

        process_stage: with a value of 'process'
        process_type: with a vlue of 'cset'
        compound_set_name: the compound set name (its primary key)

    Otherwise (i.e. not validated)...

        Returns the response from 'validate_compound_set()'
    """

    process_stage, process_type, validate_dict, validated, params = validate_output

    logger.info('process_compound_set() ENTER')
    logger.info('process_compound_set() process_type=%s', process_type)
    logger.info('process_compound_set() validated=%s', validated)
    logger.info('process_compound_set() params=%s', params)

    if not validated:
        logger.warning('process_compound_set() EXIT params=%s (not validated)', params)
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

    logger.info('process_compound_set() EXIT (CompoundSet.name="%s")', compound_set.name)
    return {'process_stage': 'process',
            'process_type': 'cset',
            'compound_set_name': compound_set.name}


@shared_task
def validate_compound_set(task_params):
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

    # Required entry params
    user_id = task_params['user_id']
    sdf_file = task_params['sdf_file']
    target = task_params['target']
    # Optional params
    zfile = task_params.get('zfile')
    update = task_params.get('update')

    logger.info('validate_compound_set() ENTER')
    logger.info('validate_compound_set() user_id=%s', user_id)
    logger.info('validate_compound_set() sdf_file=%s', sdf_file)
    logger.info('validate_compound_set() target=%s', target)
    logger.info('validate_compound_set() update=%s', update)

    validated = True
    validate_dict = {'molecule_name': [],
                     'field': [],
                     'warning_string': []}

    # params = {
    #     'sdf': 'tests/test_data/test_chodera/test_0_2.sdf',
    #     'target': 'Mpro',
    #     'choice': 1,
    #     'update_set': None,
    #     'pdb_zip': '/code/tests/test_data/test_chodera/references.zip'
    # }

    outbound_params = {
        'user_id': user_id,
        'sdf': sdf_file,
        'target': target,
        'pdb_zip': zfile
    }

    # Protect ourselves from an empty, blank or missing SD file.
    if sdf_file is None or len(sdf_file) == 0:
        validated = False
        logger.info('validate_compound_set() EXIT (no file) validated=%s outbound_params=%s',
                    validated, outbound_params)
        return 'validate', 'cset', validate_dict, validated, outbound_params
    elif not os.path.isfile(sdf_file):
        validated = False
        logger.info('validate_compound_set() EXIT (missing file) validated=%s outbound_params=%s',
                    validated, outbound_params)
        return 'validate', 'cset', validate_dict, validated, outbound_params

    suppl = Chem.SDMolSupplier(sdf_file)
    # print('%d mols detected (including blank mol)' % (len(suppl),))
    blank_mol = suppl[0]

    if blank_mol is None:
        validate_dict = add_warning(molecule_name='Blank Mol',
                                    field='N/A',
                                    warning_string='your blank molecule could not be'
                                                   ' read by rdkit. The molecule must'
                                                   ' have at least one atom!'
                                                   ' No other checks were done',
                                    validate_dict=validate_dict)
        validated = False
        logger.warning('validate_compound_set() EXIT'
                       ' user_id=%s sdf_file=%s validated=False',
                       user_id, sdf_file)
        # Can't get submitter name or method when there is now mol
        submitter_name = ''
        submitter_method = ''
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
                molecule_name = 'Unknown (no _Name property)'
                if mol.HasProp('_Name'):
                    molecule_name = mol.GetProp('_Name')
                add_warning(molecule_name=molecule_name,
                            field='property (missing)',
                            warning_string=f'{diff} property is missing from this molecule',
                            validate_dict=validate_dict)

    # Check version in blank mol
    validate_dict = check_ver_name(blank_mol, SDF_VERSION, validate_dict)

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
            molecule_name = ''
            if m.HasProp('_Name'):
                molecule_name = m.GetProp('_Name')
            validate_dict = check_name_characters(molecule_name, validate_dict)
            # validate_dict = check_pdb(m, validate_dict, target, zfile)
            validate_dict = check_refmol(m, validate_dict, target)
            validate_dict = check_field_populated(m, validate_dict)
            validate_dict = check_SMILES(m, validate_dict)

    len_validate_dict = len(validate_dict['molecule_name'])
    if len_validate_dict != 0:
        logger.warning('validate_compound_set()'
                       ' user_id=%s sdf_file=%s len_validate_dict=%s validated=False',
                       user_id, sdf_file, len_validate_dict)
        validated = False

    logger.info('validate_compound_set() EXIT validated=%s outbound_params=%s',
                validated, outbound_params)
    return 'validate', 'cset', validate_dict, validated, outbound_params


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
    logger.info('+ TASK target_zip=%s', target_zip)

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

    logger.info('- TASK target_zip=%s validated=%s validate_dict=%s',
                target_zip, validated, validate_dict)

    return ('validate', 'tset', validate_dict, validated, new_data_folder, target, proposal,
            submitter_name, email)


@celery_app.task(bind=True)
def task_load_target(self, data_bundle, proposal_ref=None, contact_email=None, user_id=None):
    logger.info('TASK load_target lauched, target_zip=%s', data_bundle)
    try:
        load_target(data_bundle, proposal_ref=proposal_ref, contact_email=contact_email, user_id=user_id, task=self)
    except KeyError as err:
        logger.error(err.args[0])
        self.update_state(
                state="ERROR",
                meta={"description": err.args[0],},
            )
    except IntegrityError as err:
        logger.error(err.args[0])
        self.update_state(
                state="ERROR",
                meta={"description": err.args[0],},
            )
    except FileNotFoundError as err:
        logger.error(err.args[0])
        self.update_state(
                state="ERROR",
                meta={"description": err.args[0],},
            )

    logger.info('TASK load_target completed, target_zip=%s', data_bundle)



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
    logger.info('+ TASK validate_output=%s', validate_output)

    # Validate output is a tuple - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    process_stage, process_type, validate_dict, validated, new_data_folder, target_name, proposal_ref, submitter_name, \
        contact_email = validate_output

    # If there is a validation error, stop here.
    if not validated:
        logger.warning('- TASK Leaving (not validated)')
        return process_stage, 'tset', validate_dict, validated

    logger.info('+ TASK Calling process_target(%s, %s, %s)...',
                new_data_folder, target_name, proposal_ref)
    mols_loaded, mols_processed = process_target(new_data_folder, target_name, proposal_ref)
    logger.info('+ TASK process_target() returned mols_loaded=%s mols_processed=%s',
                mols_loaded, mols_processed)

    return 'process', 'tset', target_name, mols_loaded, mols_processed, contact_email


@shared_task
def process_job_file_transfer(auth_token, jt_id):
    """ Celery task to take a list of proteins and specification and transfer the files to Squonk2

    Parameters
    ----------
    task_
        jt_id of job_file_transfer record

    Returns
    -------
    target str
        name of the processed target set

    """

    logger.info('+ TASK Starting File Transfer (%s) [STARTED]', id)
    job_transfer = JobFileTransfer.objects.get(id=id)
    job_transfer.transfer_status = "STARTED"
    job_transfer.transfer_task_id = str(process_job_file_transfer.request.id)
    job_transfer.save()
    try:
        process_file_transfer(auth_token, job_transfer.id)
    except RuntimeError as error:
        logger.error('- TASK File transfer (%s) [FAILED] error=%s',
                     id, error)
        logger.error(error)
        job_transfer.transfer_status = "FAILURE"
        job_transfer.save()
    else:
        # Update the transfer datetime for comparison with the target upload datetime.
        # This should only be done on a successful upload.
        job_transfer.transfer_datetime = datetime.datetime.now(datetime.timezone.utc)
        job_transfer.transfer_progress = 100.00
        job_transfer.transfer_status = "SUCCESS"
        files_spec = {"proteins": list(SQUONK_PROT_MAPPING.keys()),
                      "compounds": list(SQUONK_COMP_MAPPING.keys())}
        job_transfer.transfer_spec = files_spec
        job_transfer.save()
        logger.info('+ TASK File transfer (%s) [SUCCESS]', id)

    return job_transfer.transfer_status


@shared_task
def process_compound_set_job_file(task_params):
    """Celery task to retrieve files generated by a JobRequest on Squonk2.
    This is invoked if the job definition supports/requires uploading of Job results.

    Parameters
    ----------
    task_params:
        jr_id: id of job_request record
        transition_time: The time the JOb ran

    Returns
    -------
    parameters suitable for validation
    """

    jr_id = task_params['jr_id']
    transition_time = task_params['transition_time']
    job_output_path = task_params['job_output_path']
    job_output_filename = task_params['job_output_filename']

    logger.info('+ TASK task_params=%s', task_params)

    job_request = JobRequest.objects.get(id=jr_id)
    job_request.upload_task_id = str(process_compound_set_job_file.request.id)
    job_request.save()

    sd_file = None
    try:
        sd_file = process_compound_set_file(jr_id,
                                            transition_time,
                                            job_output_path,
                                            job_output_filename)
    except RuntimeError as error:
        logger.error('- TASK file Upload failed (%s)', jr_id)
        logger.error(error)
    else:
        # Update the transfer datetime for comparison with the target upload datetime.
        # This should only be done on a successful upload.
        logger.info('+ TASK file Upload Ended Successfully (%s)', jr_id)

    # We are expected to be followed by 'validate_compound_set'
    # which expects user_id, sdf_file and target
    return {'user_id': job_request.user.id,
            'sdf_file': sd_file,
            'target': job_request.target.title}


@shared_task
def erase_compound_set_job_material(task_params, job_request_id=0):
    """Celery task to clean-up files generated by a JobRequest on Squonk2.
    We receive the output of 'process_compound_set()'. If the first field is not
    'process' then we can assume the upload failed, maybe during validation?

    Parameters
    ----------
    task_params:
        Consists of result of process_compound_set()

    Returns
    -------
    parameters suitable for validation
    """
    assert task_params
    assert len(task_params) > 0

    # Do nothing if the job-request isn't set.
    if not job_request_id:
        return

    job_request = JobRequest.objects.get(id=job_request_id)
    logger.info('+ TASK Erasing job material job_request %s', job_request)

    # Upload done (successfully or not)
    # 'task_params' (a dictionary) is expected to be
    # the return value of 'process_compound_set()'
    # so set the upload status. We expect to find
    # 'process_stage' and 'compound_set_name'.
    #
    # Task linking is a bit of a mess atm,
    # if something went wrong we'll get a tuple, not a dictionary.
    if isinstance(task_params, dict) \
            and task_params['process_stage'] == 'process' \
            and task_params['compound_set_name']:
        logger.info('Upload successful (%d) CompoundSet.name="%s"',
                    job_request_id, task_params['compound_set_name'])
        job_request.upload_status = 'SUCCESS'
        # We're given a compound set name.
        # Get its record and put that into the JobRequest...
        cs = ComputedSet.objects.get(name=task_params['compound_set_name'])
        assert cs
        job_request.computed_set = cs
    else:
        # Failed validation?
        logger.info('- TASK Upload failed (%d) - task_params=%s',
                    job_request_id, task_params)
        logger.warning('Upload failed (%d) - process_stage value not satisfied',
                       job_request_id)
        job_request.upload_status = 'FAILURE'
    job_request.save()
    logger.info('+ TASK Erased and updated job_request %s', job_request)

    # Always erase uploaded data
    delete_media_sub_directory(get_upload_sub_directory(job_request))
