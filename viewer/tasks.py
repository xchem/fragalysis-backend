import os
import datetime
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
django.setup()
from django.conf import settings

from rdkit import Chem
from rdkit.Chem import Descriptors
from viewer.models import Target, Compound, DesignSet
import os.path
import zipfile

from celery import shared_task
from .sdf_check import *
from .compound_set_upload import *
from .target_set_upload import process_target, validate_target

# import the logging library
#import logging
# Get an instance of a logger
#logger = logging.getLogger(__name__)

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
def add_cset_mols(cset, sdf_file, target=None, zfile=None):
    # check the computed set exists
    # for each molecule:
    # - check there's a ScoreDescription for each score
    # - add the molecule
    target_obj = Target.objects.get(title=target)

    existing = ComputedSet.objects.filter(unique_name=cset, target=target_obj)

    if len(existing) == 1:
        compound_set = existing[0]
    else:
        return 'Computed set does not exist, check its name or upload it at viewer/upload_cset'

    mols = get_additional_mols(filename=sdf_file, compound_set=compound_set)

    # if mols is a string (error) return it
    if isinstance(mols, str):
        return mols

    # process every other mol
    for i in range(0, len(mols)):
        process_mol(mol=mols[i], target=target, compound_set=compound_set, filename=sdf_file, zfile=zfile)

    # need to add new mols to old sdf!
    old_sdf = str(compound_set.submitted_sdf.path)
    with open(old_sdf, 'ab') as w:
        w.write(b'\n')
        w.write(open(sdf_file, 'rb').read())

    return f"{compound_set.unique_name} was successfully updated"


@shared_task
def process_compound_set(validate_output):
    """ Celery task to process a computed set, that takes the output of the validation task, and uploads molecules to a
    new computed set if the uploaded files are valid

    Parameters
    ----------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded sdf file
            - target (str): name of the target that the computed set is associated with
            - zfile (dict): dictionary where key is the name of the file minus extension and path, and value is the filename, which is saved to temporary storage by `viewer.views.UploadCSet`
            - submitter_name (str): name of the author of the computed set
            - submitter_method (str): name of the method used to generate the computed set

    Returns
    -------
    compound_set.name: str
        name of the computed set

    """
    # Validate output is a tuple - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    process_type, validate_dict, validated, filename, target, zfile, \
    submitter_name,  submitter_method = validate_output

    if not validated:
        return (validate_dict, validated)

    if validated:
        #print('processing compound set: ' + filename)
        filename = str(filename)

        # create a new compound set
        set_name = ''.join(filename.split('/')[-1].replace('.sdf','').split('_')[1:])

        existing = ComputedSet.objects.filter(unique_name="".join(submitter_name.split()) + '-' + "".join(submitter_method.split()))

        if len(existing)==1:
            compound_set=existing[0]
        else:
            compound_set = ComputedSet()

        compound_set.name = set_name
        matching_target = Target.objects.get(title=target)
        compound_set.target = matching_target
        ver = float(version.strip('ver_'))
        compound_set.spec_version = ver
        compound_set.unique_name = "".join(submitter_name.split()) + '-' + "".join(submitter_method.split())
        compound_set.save()

        # set descriptions and get all other mols back
        mols_to_process = set_descriptions(filename=filename, compound_set=compound_set)

        # process every other mol
        for i in range(0, len(mols_to_process)):
            process_mol(mols_to_process[i], target, compound_set, filename, zfile)

        # check that molecules have been added to the compound set
        check = ComputedMolecule.objects.filter(computed_set=compound_set)
        #print(str(len(check)) + '/' + str(len(mols_to_process)) + ' succesfully processed in ' + set_name + ' cpd set')

        # move and save the compound set
        new_filename = settings.MEDIA_ROOT + 'compound_sets/' + filename.split('/')[-1]
        os.rename(filename, new_filename)
        compound_set.submitted_sdf = new_filename
        compound_set.save()

        # if no molecules were processed, delete the compound set
        if len(check) == 0:
            compound_set.delete()
            #print('No molecules processed... deleting ' + set_name + ' compound set')
            return None

        return 'process', 'cset', compound_set.name

# End Uploading Compound Sets ###

# Validating Compound Sets ###

# Set .sdf format version here
version = 'ver_1.2'


@shared_task
def validate_compound_set(sdf_file, target=None, zfile=None, update=None):
    """ Celery task to process validate the uploaded files for a computed set upload. SDF file is mandatory, zip file is
    optional

    Parameters
    ----------
    sdf_file: str
        filepath of the uploaded sdf file, which is saved to temporary storage by `viewer.views.UploadCSet`
    target: str
        name of the target (`viewer.models.Target.title`) to add add the computed set to
    zfile: dict
        dictionary where key is the name of the file minus extension and path, and value is the filename, which is saved
        to temporary storage by `viewer.views.UploadCSet`

    Returns
    -------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the calidation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded sdf file
            - target (str): name of the target that the computed set is associated with
            - zfile (dict): dictionary where key is the name of the file minus extension and path, and value is the
              filename, which is saved to temporary storage by `viewer.views.UploadCSet`
            - submitter_name (str): name of the author of the computed set
            - submitter_method (str): name of the method used to generate the computed set

    """
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
                                    warning_string='your blank molecule could not be read by rdkit. The molecule must have at least one atom! No other checks were done',
                                    validate_dict=validate_dict)
        validated = False
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
                        warning_string='SDF entry number: %s can not be converted into an rdkit mol object' % (index,),
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
                            warning_string='%s property is missing from this molecule' % (diff,),
                            validate_dict=validate_dict)

    # Check version in blank mol
    validate_dict = check_ver_name(blank_mol, version, validate_dict)

    # Check compulsory fields in blank mol props
    validate_dict = check_blank_mol_props(blank_mol, validate_dict)

    # Check properties have been described and validate url
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
            validate_dict = check_pdb(m, validate_dict, target, zfile)
            validate_dict = check_refmol(m, validate_dict, target)
            validate_dict = check_field_populated(m, validate_dict)
            validate_dict = check_SMILES(m, validate_dict)

    if len(validate_dict['molecule_name']) != 0:
        validated = False

    unique_name = "".join(submitter_name.split()) + '-' + "".join(submitter_method.split())
    csets = ComputedSet.objects.filter(unique_name=unique_name)
    [c.delete() for c in csets]

    return ('validate', validate_dict, validated, sdf_file, target, zfile,
            submitter_name,  submitter_method)

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
        molecules = Molecule.objects.filter(prot_id__code__contains=insp.split('_')[0])
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
    for i, row in set_df.iterrows():
        compounds.append(process_design_compound(row))

    for compound in compounds:
        new_set.compounds.add(compound)

    new_set.save()

    return compounds


def process_design_sets(df, set_type=None, set_description=None):
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
def validate_target_set(target_zip, target=None, update=None, proposal=None):
    """ Celery task to process validate the uploaded files/format for a target set upload. Zip file is mandatory

    Parameters
    ----------
    target_zip: str
        filepath of the uploaded target file, which is saved to temporary storage by `viewer.views.UploadTSet`
    target: str
        name of the target (`viewer.models.Target.title`) to add add the target set to
    update: dict
        dictionary where key is the name of the file minus extension and path, and value is the filename, which is
        saved to temporary storage by `viewer.views.UploadTSet`

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
    # This will create the target folder in the tmp/ location.
    with zipfile.ZipFile(target_zip, 'r') as zip_ref:
        zip_ref.extractall(tmp_folder)

    validated, validate_dict = validate_target(tmp_folder, target, proposal)

    os.remove(target_zip)

    return ('validate', validate_dict, validated, tmp_folder, target, proposal,
            submitter_name)


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
    process_type, validate_dict, validated, tmp_folder, target_name, proposal_ref, submitter_name = validate_output

    # If there is a validation error, stop here.
    if not validated:
        return process_type, validate_dict, validated

    if validated:
        logger.info('+ processing target set: ' + target_name + ' target_folder:' + tmp_folder)
        process_target(tmp_folder, target_name, proposal_ref)
        return 'process', 'tset', target_name

# End Target Sets ###
