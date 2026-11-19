import collections
import datetime
import logging
import os
import shutil

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")

import django

django.setup()

import zipfile

from celery import shared_task
from celery.utils.log import get_task_logger
from django.conf import settings
from rdkit import Chem
from rdkit.Chem import Descriptors

from fragalysis.celery import app as celery_app
from viewer.models import Compound, DesignSet
from viewer.target_loader import load_target

from .cache import clear_view_cache
from .cset_upload import MolOps, PdbOps, blank_mol_vals
from .download_structures import create_download
from .models import ComputedSet, JobFileTransfer, JobRequest, SiteObservation
from .squonk_job_file_transfer import process_file_transfer
from .squonk_job_file_upload import get_upload_sub_directory, process_compound_set_file
from .target_set_upload import validate_target
from .utils import delete_media_sub_directory

# If Celery configured to always run 'synchronously' (eager),
# then use a standard logger, otherwise use the Celery logger.
if settings.CELERY_TASK_ALWAYS_EAGER:
    logger = logging.getLogger(__name__)
else:
    logger = get_task_logger(__name__)


@shared_task
def process_compound_set(task_params):
    """Celery task to upload a computed set. Started by views.py's
    UploadComputedSetView (and the Squonk job chain), it loads the molecules
    into a new (or updated) computed set. SDF file is mandatory, zip file
    is optional.

    There is no separate validation step: anything wrong with the files
    raises, the upload is rolled back and the task fails with the error.

    (Code is legacy, and has been taken on by multiple users - comments often
    retro-fitted and no apologies for the design)

    Parameters
    ----------
    task_params: dict
        - user_id (int): User record ID of user initiating the task (may be None)
        - sdf_file (str): path of the uploaded sdf file, saved to temporary
          storage by `viewer.views.UploadComputedSetView`
        - target (int): id of the target the computed set is associated with
        - zfile (str): optional, path of the uploaded zip of PDB files
        - update (str): optional, the name of a computed set to update, or 'None'

    Returns
    -------
    A tuple that contains the fields: -

        A value of 'process'
        the compound set name (its primary key)
        a dict of any warnings collected while processing
    """

    user_id = task_params['user_id']
    sdf_file = task_params['sdf_file']
    target = task_params['target']
    pdb_zip = task_params.get('zfile')
    computed_set_id = task_params.get('update')

    logger.info('process_compound_set() ENTER')
    logger.info('process_compound_set() task_params=%s', task_params)

    # Protect ourselves from an empty, blank or missing SD file.
    if not sdf_file or not os.path.isfile(sdf_file):
        msg = f'SD file "{sdf_file}" is missing'
        logger.warning('process_compound_set() EXIT (%s)', msg)
        raise ValueError(msg)

    submitter_name, submitter_method, blank_version = blank_mol_vals(sdf_file)
    zfile, zfile_hashvals = PdbOps().save_pdb_zip(pdb_zip)

    # We expect a User ID, but it may be None.
    # If AUTHENTICATE_UPLOAD is False the User ID will be None.
    save_mols = MolOps(
        user_id=user_id,
        sdf_filename=sdf_file,
        submitter_name=submitter_name,
        submitter_method=submitter_method,
        target=target,
        version=blank_version,
        zfile=zfile,
        zfile_hashvals=zfile_hashvals,
        computed_set_id=computed_set_id,
    )
    compound_set, process_messages = save_mols.task()

    # ComputedMolecule rows have just been written; drop both
    # ComputedMoleculesView and ComputedMolAndScoreView caches.
    clear_view_cache("computed-molecules")

    logger.info('process_compound_set() EXIT (CompoundSet.id="%s")', compound_set.id)
    return 'process', compound_set.id, process_messages


def create_mol(inchi, long_inchi=None, name=None):
    # check for an existing compound
    if long_inchi:
        cpd = Compound.objects.filter(long_inchi=long_inchi)
        sanitized_mol = Chem.MolFromInchi(long_inchi, sanitize=True)
    else:
        cpd = Compound.objects.filter(inchi=inchi)
        sanitized_mol = Chem.MolFromInchi(inchi, sanitize=True)

    new_mol = cpd[0] if len(cpd) != 0 else Compound()
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
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    sanitized_mol_smiles = Chem.MolToSmiles(mol, canonical=True)
    sanitized_mol = Chem.MolFromSmiles(sanitized_mol_smiles)
    inchi = Chem.inchi.MolToInchi(sanitized_mol)

    if len(inchi) > 255:
        # TODO: get_inchi in model
        long_inchi = inchi
        inchi = str(long_inchi)[:255]
        logger.warning(
            'inchi too long, truncating to 255 characters (original inchi=%s, sanitized_mol=%s, new inchi=%s)',
            long_inchi,
            sanitized_mol,
            inchi,
        )

    new_mol = create_mol(inchi)

    # deal with inspirations
    inspirations = compound_row['inspirations'].split(',')
    for insp in inspirations:
        # TODO: find matching molecules - change to molecules and search history to find the correct version.
        #  -- search all history and find most recent with matching code? or code most closely matching design date?
        # (Can this be accessed, or does the view need changing to display the correct one? Not implemented yet anyway)
        # molecules = Molecule.objects.filter(prot_id__code__contains=insp.split(':')[0].split('_')[0])
        molecules = SiteObservation.objects.filter(
            prot_id__code__contains=insp.split(':')[0].split('_')[0]
        )
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

    compounds = [process_design_compound(row) for _, row in set_df.iterrows()]
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
    """Celery task to process validate the uploaded files/format for a target set upload. Zip file is mandatory

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

    logger.info(
        '- TASK target_zip=%s validated=%s validate_dict=%s',
        target_zip,
        validated,
        validate_dict,
    )

    return (
        'validate',
        'tset',
        validate_dict,
        validated,
        new_data_folder,
        target,
        proposal,
        submitter_name,
        email,
    )


@celery_app.task(bind=True)
def task_load_target(
    self, data_bundle=None, proposal_ref=None, user_id=None, curation_file=None
):
    logger.info(
        'TASK %s load_target launched, target_zip=%s', self.request.id, data_bundle
    )
    load_target(
        data_bundle,
        proposal_ref=proposal_ref,
        user_id=user_id,
        task=self,
        curation_file=curation_file,
    )
    logger.info(
        'TASK %s load_target completed, target_zip=%s', self.request.id, data_bundle
    )


@shared_task
def process_job_file_transfer(auth_token, jt_id):
    """Celery task to take a list of proteins and a specification
    and transfer the files to Squonk2

    Parameters
    ----------
    auth_token for the executing user
    jt_id is the job_file_transfer record ID

    Returns
    -------
    target str
        name of the processed target set

    """

    logger.info('+ File Transfer (id=%s) [STARTED]', jt_id)
    job_transfer = JobFileTransfer.objects.get(id=jt_id)
    job_transfer.transfer_status = "STARTED"
    job_transfer.transfer_task_id = str(process_job_file_transfer.request.id)
    job_transfer.save()

    try:
        process_file_transfer(auth_token, job_transfer.id)
    except RuntimeError as error:
        job_transfer.transfer_status = "FAILURE"
        job_transfer.transfer_datetime = datetime.datetime.now(datetime.timezone.utc)
        job_transfer.save()
        logger.error('- File transfer (id=%s) [RuntimeError "%s"]', jt_id, error)
    else:
        # Update the transfer datetime for comparison with the target upload datetime.
        # This should only be done on a successful upload.
        job_transfer.transfer_datetime = datetime.datetime.now(datetime.timezone.utc)
        job_transfer.transfer_progress = 100.00
        job_transfer.transfer_status = "SUCCESS"
        job_transfer.save()
        logger.info('- File transfer (id=%s) [SUCCESS]', jt_id)

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
    parameters suitable for process_compound_set()
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
        sd_file = process_compound_set_file(
            jr_id, transition_time, job_output_path, job_output_filename
        )
    except RuntimeError as error:
        logger.error('- TASK file Upload failed (%s)', jr_id)
        logger.error(error)
    else:
        # Update the transfer datetime for comparison with the target upload datetime.
        # This should only be done on a successful upload.
        logger.info('+ TASK file Upload Ended Successfully (%s)', jr_id)

    # We are expected to be followed by 'process_compound_set'
    # which expects user_id, sdf_file and target
    return {
        'user_id': job_request.user.id,
        'sdf_file': sd_file,
        'target': job_request.target.pk,
    }


@shared_task
def erase_compound_set_job_material(task_params, job_request_id=0):
    """Celery task to clean-up files generated by a JobRequest on Squonk2.
    We receive the output of 'process_compound_set()'. It provides us with a tuple.
    If the first field is not 'process' then we can assume the upload failed.

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

    # Upload done (successfully or not).
    # 'task_params' (a list/tuple) is expected to be the return value of
    # 'process_compound_set()' so set the upload status. We expect to find
    # 'process_stage' [index 0] and 'compound_set.id' [index 1].
    #

    # Inter-task parameter linking is a bit of a mess atm,
    # for now deal with tuple, or list (i.e. any "Sequence").
    #
    # Assume failure...
    job_request.upload_status = 'FAILURE'
    if (
        isinstance(task_params, collections.abc.Sequence)
        and not isinstance(task_params, str)
        and len(task_params) > 1
    ):
        process_stage = task_params[0]
        if process_stage == 'process':
            # We've come from the right task...
            cs_id: str = task_params[1]
            logger.info(
                'Upload successful (%d) ComputedSet.id="%s"', job_request_id, cs_id
            )
            job_request.upload_status = 'SUCCESS'
            # We're given a compound set is.
            # Get its record and put that into the JobRequest...
            cs = ComputedSet.objects.get(pk=cs_id)
            assert cs
            job_request.computed_set = cs
        else:
            logger.warning(
                "Upload failed (%d) - incorrect process_stage. Expected 'process', got '%s'."
                " task_params=%s",
                job_request_id,
                process_stage,
                task_params,
            )
    else:
        # Invalid task_params
        logger.info(
            "- TASK Upload failed (%d) - unexpected task_params value."
            " Need a Sequence of more than 1 item, got '%s'",
            job_request_id,
            task_params,
        )

    assert job_request.upload_status in {'SUCCESS', 'FAILURE'}
    job_request.save()
    logger.info('+ TASK Erased and updated JobRequest %s', job_request)

    # Always erase uploaded data
    delete_media_sub_directory(get_upload_sub_directory(job_request))


@celery_app.task(bind=True)
def task_create_download(
    self,
    download_link_id: int,
    use_zip: bool = False,
):
    logger.info(
        'TASK %s create_download_link launched, download_link=%s',
        self.request.id,
        download_link_id,
    )
    create_download(download_link_id, task=self, use_zip=use_zip)
    logger.info('TASK %s create_download_link completed', self.request.id)
