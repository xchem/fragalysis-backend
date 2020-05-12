import os
import datetime
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
django.setup()
from django.conf import settings

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from rdkit import Chem
from viewer.models import (
    Target,
    CompoundSet,
    ComputedCompound,
    ScoreDescription,
    NumericalScoreValues,
    TextScoreValues,
    Protein,
    Molecule,
    CompoundSetSubmitter)
import ast
import os.path


import validators
import numpy as np
import os

from celery import shared_task
import psutil

### Uploading ###

def dataType(str):
    str = str.strip()
    if len(str) == 0: return 'BLANK'
    try:
        t = ast.literal_eval(str)

    except ValueError:
        return 'TEXT'
    except SyntaxError:
        return 'TEXT'

    else:
        if type(t) in [int, long, float, bool]:
            if t in {True, False, 'TRUE', 'FALSE', 'true', 'false', 'yes', 'no', 'YES', 'NO', 'Yes', 'No'}:
                return 'BIT'
            if type(t) is int or type(t) is long:
                return 'INT'
            if type(t) is float:
                return 'FLOAT'
        else:
            return 'TEXT'

def get_inspiration_frags(cpd, compound_set):
    pass

# use zfile object for pdb files uploaded in zip
def get_prot(mol, compound_set, zfile):
    pdb_option = mol.GetProp('ref_pdb')
    name = pdb_option.split('/')[-1]
    if zfile:
        if pdb_option in zfile['zf_list']:
            data = zfile['zip_obj'].read(pdb_option)
            field = default_storage.save('tmp/' + name, ContentFile(data))

    else:
        name = compound_set.target.title + '-' + pdb_option
        print('PROT: ' + name)
        prot = Protein.objects.get(code__contains=name)
        field = prot.pdb_info

    return field


def set_props(cpd, props, compound_set):
    if 'ref_mols' and 'ref_pdb' not in props.keys():
        raise Exception('ref_mols and ref_pdb not set!')
    set_obj = ScoreDescription.objects.filter(compound_set=compound_set)
    set_props_list = [s.name for s in set_obj]
    for key in props.keys():
        if key in set_props_list not in ['ref_mols', 'ref_pdb', 'original SMILES']:
            if dataType(str(props[key]))=='TEXT':
                score_value = TextScoreValues()
            else:
                score_value = NumericalScoreValues()
            score_value.score = ScoreDescription.objects.get(compound_set=compound_set,
                                                             name=key)
            score_value.value = props[key]
            score_value.compound = cpd
            score_value.save()

    return set_obj


def set_mol(mol, compound_set, filename, zfile=None):
    # zfile = {'zip_obj': zf, 'zf_list': zip_names}

    smiles = Chem.MolToSmiles(mol)
    name = mol.GetProp('_Name')
    mol_block = Chem.MolToMolBlock(mol)

    insp = mol.GetProp('ref_mols')
    insp = insp.split(',')
    insp = [i.strip() for i in insp]
    insp_frags = []
    for i in insp:
        mols = Molecule.objects.filter(prot_id__code__contains=str(compound_set.target.title + '-' + i),
                                       prot_id__target_id=compound_set.target)
        if len(mols)>1:
            ids = [m.cmpd_id.id for m in mols]
            ind = ids.index(max(ids))
            ref = mols[ind]
        if len(mols)==1:
            ref = mols[0]
        if len(mols)==0:
            raise Exception('No matching molecules found for inspiration frag ' + i)

        insp_frags.append(ref)


    orig = mol.GetProp('original SMILES')

    prot_field = get_prot(mol, compound_set, zfile)
    if 'tmp' in prot_field:
        # move and save the compound set
        old_filename = settings.MEDIA_ROOT + prot_field
        new_filename = settings.MEDIA_ROOT + 'pdbs/' + prot_field.split('/')[-1]
        os.rename(old_filename, new_filename)
        prot_field = new_filename
        # compound_set.save()

    cpd = ComputedCompound()
    cpd.sdf_info = mol_block
    cpd.compound_set = compound_set
    cpd.name = name
    cpd.smiles = smiles
    cpd.pdb_info = prot_field
    cpd.original_smiles = orig

    cpd.save()

    [cpd.inspiration_frags.add(mol) for mol in insp_frags]

    cpd.save()

    return cpd


def process_mol(mol, compound_set, filename, zfile=None):
    cpd = set_mol(mol, compound_set, filename, zfile)
    other_props = mol.GetPropsAsDict()
    compound_set = set_props(cpd, other_props, compound_set)

    return compound_set


def get_submission_info(description_mol):
    y_m_d = description_mol.GetProp('generation_date').split('-')

    submitter_dict = {'name': description_mol.GetProp('submitter_name'),
                      'email': description_mol.GetProp('submitter_email'),
                      'institution': description_mol.GetProp('submitter_institution'),
                      'generation_date': datetime.date(int(y_m_d[0]), int(y_m_d[1]), int(y_m_d[2])),
                      'method': description_mol.GetProp('method')}

    submitter = CompoundSetSubmitter(**submitter_dict)
    submitter.save()
    return submitter


def set_descriptions(filename, compound_set):
    suppl = Chem.SDMolSupplier(str(filename))
    description_mol = suppl[0]

    submitter = get_submission_info(description_mol)

    description_dict = description_mol.GetPropsAsDict()
    version = description_mol.GetProp('_Name')
    compound_set.spec_version = version.split('_')[-1]
    method = description_mol.GetProp('ref_url')
    compound_set.method_url = method
    compound_set.submitter = submitter
    compound_set.save()

    for key in description_dict.keys():
        desc = ScoreDescription()
        desc.compound_set = compound_set
        desc.name = key
        desc.description = description_dict[key]
        desc.save()

    mols = []
    for i in range(1, len(suppl)):
        mols.append(suppl[i])

    return mols

def check_services():
    services = [p.name() for p in psutil.process_iter()]
    if 'redis-server' not in services:
        os.system('redis-server &')
    if 'celery' not in services:
        os.system('celery -A fragalysis worker -l info &')
    services = [p.name() for p in psutil.process_iter()]
    if 'redis-server' not in services or 'celery' not in services:
        return False
    return True


@shared_task(bind=True)
def process_compound_set(self, target, filename, zfile=None):
    print('processing compound set: ' + filename)
    filename = str(filename)
    # create a new compound set
    set_name = ''.join(filename.split('/')[-1].replace('.sdf','').split('_')[1:])
    compound_set = CompoundSet()
    compound_set.name = set_name
    matching_target = Target.objects.get(title=target)
    compound_set.target = matching_target

    # set descriptions and get all other mols back
    mols_to_process = set_descriptions(filename=filename, compound_set=compound_set)

    # process every other mol
    for i in range(0, len(mols_to_process)):
        process_mol(mols_to_process[i], compound_set, filename, zfile)

    # check that molecules have been added to the compound set
    check = ComputedCompound.objects.filter(compound_set=compound_set)
    print(str(len(check)) + '/' + str(len(mols_to_process)) + ' succesfully processed in ' + set_name + ' cpd set')

    # move and save the compound set
    new_filename = settings.MEDIA_ROOT + 'compound_sets/' + filename.split('/')[-1]
    os.rename(filename, new_filename)
    compound_set.submitted_sdf = new_filename
    compound_set.save()

    # if no molecules were processed, delete the compound set
    if len(check) == 0:
        compound_set.delete()
        print('No molecules processed... deleting ' + set_name + ' compound set')
        return None

    return compound_set.name

### End Uploading ###

### Validating ###

# Set .sdf format version here
version = 'ver_1.2'


def check_compound_set(description_mol, validate_dict):
    y_m_d = description_mol.GetProp('generation_date').split('-')

    submitter_dict = {'submitter__name': description_mol.GetProp('submitter_name'),
                      'submitter__email': description_mol.GetProp('submitter_email'),
                      'submitter__institution': description_mol.GetProp('submitter_institution'),
                      'submitter__generation_date': datetime.date(int(y_m_d[0]), int(y_m_d[1]), int(y_m_d[2])),
                      'submitter__method': description_mol.GetProp('method')}

    query = CompoundSet.objects.filter(**submitter_dict)

    if len(query) != 0:
        validate_dict = add_warning(molecule_name='File error',
                                    field='compound set',
                                    warning_string="a compound set with the auto_generated name " + query[
                                        0].submitter.unique_name + " already exists (change method name in blank mol method field)",
                                    validate_dict=validate_dict)

    return validate_dict


def add_warning(molecule_name, field, warning_string, validate_dict):
    validate_dict['molecule_name'].append(molecule_name)
    validate_dict['field'].append(field)
    validate_dict['warning_string'].append(warning_string)

    return validate_dict


def check_sdf(sdf_file, validate_dict):
    """
    Checks if .sdf file can be read and follows naming format:
    'compound-set_<name>.sdf' with <name> replaced with
    the name you wish to give it. e.g. compound-set_fragmenstein.sdf

    :sdf_file: is the sdf in the specified format
    :return: Updates validate dictionary with pass/fail message
    """
    # Check filename
    if sdf_file.startswith("compound-set_") and sdf_file.endswith(".sdf") is False:
        validate_dict = add_warning(molecule_name='File error',
                                    field='_File_name',
                                    warning_string="illegal filename: " + str(sdf_file) + " found",
                                    validate_dict=validate_dict)

    return validate_dict


def check_refmol(mol, validate_dict, target=None):
    if target:
        refmols = mol.GetProp('ref_mols').split(',')
        for ref in refmols:
            query = Protein.objects.filter(code__contains=target + '-' + ref.strip())
            if len(query) == 0:
                validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                            field='ref_mol',
                                            warning_string="molecule for " + str(
                                                ref.strip()) + " does not exist in fragalysis (make sure the code is exactly as it appears in fragalysis - e.g. x0123_0)",
                                            validate_dict=validate_dict)
    return validate_dict


def check_pdb(mol, validate_dict, target=None, zfile=None):
    """
    Checks if .pdb file can be read

    :mol: rdkit mol read from SD file
    :return: Updates validate dictionary with pass/fail message
    """

    # Check if pdb filepath given and exists
    test_fp = mol.GetProp('ref_pdb')

    # {'zip_obj': zf, 'zf_list': zip_names}

    if zfile:
        pdb_option = mol.GetProp('ref_pdb')
        # name = pdb_option.split('/')[-1]
        if zfile:
            if pdb_option not in zfile['zf_list']:
                validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                            field='ref_pdb',
                                            warning_string="path " + str(
                                                pdb_option) + " can't be found in uploaded zip file (list: " + str(
                                                zfile['zf_list']) + ")",
                                            validate_dict=validate_dict)

    # else:
    if target and not test_fp.endswith(".pdb"):
        query = Protein.objects.filter(code__contains=str(target + '-' + test_fp))
        if len(query) == 0:
            validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                        field='ref_pdb',
                                        warning_string="pdb for " + str(test_fp) + " does not exist in fragalysis",
                                        validate_dict=validate_dict)

    return validate_dict


def check_SMILES(mol, validate_dict):
    """
    Checks if SMILES can be read by rdkit

    :mol: rdkit mol read from SD file
    :return: Updates validate dictionary with pass/fail message
    """
    # Check SMILES
    smi_check = mol.GetProp('original SMILES')

    m = Chem.MolFromSmiles(smi_check, sanitize=False)
    if m is None:
        validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                    field='original SMILES',
                                    warning_string="invalid SMILES %s" % (smi_check,),
                                    validate_dict=validate_dict)

    return validate_dict


def check_ver_name(blank_mol, version, validate_dict):
    """
    Checks if blank mol:
    The name (title line) of this molecule should be the
    file format specification version e.g. ver_1.0 (as defined in this document)

    :blank_mol: rdkit mol of blank mol from an SD file
    :return: Updates validate dictionary with pass/fail message
    """

    ver_name = blank_mol.GetProp('_Name')
    if ver_name != version:
        validate_dict = add_warning(molecule_name=blank_mol.GetProp('_Name'),
                                    field='_Name',
                                    warning_string='illegal version: %s found. Should be %s' % (ver_name, version),
                                    validate_dict=validate_dict)

    return validate_dict


def check_blank_mol_props(mol, validate_dict):
    # check for compulsory fields in blank mols
    fields = ['ref_url', 'submitter_name', 'submitter_email', 'submitter_institution', 'generation_date', 'method']
    for field in fields:
        validate_dict = missing_field_check(mol, field, validate_dict)

    return validate_dict


def check_blank_prop(blank_mol, validate_dict):
    """
    Checks if blank mol properties have a description

    :blank_mol: rdkit mol of blank mol from an SD file
    :return: Updates validate dictionary with pass/fail message
    """

    # Check if properties populated
    property_dict = blank_mol.GetPropsAsDict()

    # Properties to ignore
    prop_ignore_list = ['ref_mols', 'ref_pdb']

    for key, value in zip(property_dict.keys(), property_dict.values()):
        if value == '' and key not in prop_ignore_list:
            validate_dict = add_warning(molecule_name=blank_mol.GetProp('_Name'),
                                        field=key,
                                        warning_string='description for %s missing' % (key,),
                                        validate_dict=validate_dict)
        if key == 'ref_url' and check_url(value) == False:
            validate_dict = add_warning(molecule_name=blank_mol.GetProp('_Name'),
                                        field=key,
                                        warning_string='illegal URL %s provided' % (value,),
                                        validate_dict=validate_dict)

    return validate_dict


def check_field_populated(mol, validate_dict):
    """
    Checks if all compulsory fields are populated:
        1. ref_mols - a comma separated list of the fragments
        2. ref_pdb - either (a) a filepath (relative to the sdf file)
            to an uploaded pdb file
        3. original SMILES - the original smiles of the compound
            before any computation was carried out

    :mol: rdkit mol other than blank_mol
    :return: Updates validate dictionary with pass/fail message
    """

    # Compuslory fields
    compulsory_fields = ['ref_pdb', 'ref_mols', 'original SMILES']

    property_dict = mol.GetPropsAsDict()
    for key, value in zip(property_dict.keys(), property_dict.values()):
        if value == '' and key in compulsory_fields:
            validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                        field=key,
                                        warning_string='value for %s missing' % (key,),
                                        validate_dict=validate_dict)

    return validate_dict


def check_url(value):
    """
    Checks if url provided exists. No internet connection required.
    Checks URL using Validators package

    :value: value associated with 'ref_url' key
    :return: False if URL can not be validated
    """

    valid = validators.url(value)
    if valid != True:
        return False


def check_name_characters(name, validate_dict):
    legal_non_alnum = ['-', '_', '.']
    for char in name:
        if not char.isalnum() and char not in legal_non_alnum:
            validate_dict = add_warning(molecule_name=name,
                                        field='_Name',
                                        warning_string='illegal character %s found' % (char,),
                                        validate_dict=validate_dict)

    return validate_dict


def missing_field_check(mol, field, validate_dict):
    props_dict = mol.GetPropsAsDict()
    if not field in props_dict.keys():
        validate_dict = add_warning(molecule_name=mol.GetProp('_Name'),
                                    field=field,
                                    warning_string='%s field not found!' % (field,),
                                    validate_dict=validate_dict)

    return validate_dict


def check_mol_props(mol, validate_dict):
    # check for missing fields
    fields = ['ref_pdb', 'ref_mols', 'original SMILES']
    for field in fields:
        validate_dict = missing_field_check(mol, field, validate_dict)

    return validate_dict

@shared_task(bind=True)
def validate(self, sdf_file, target=None, zfile=None):
    validated = True
    validate_dict = {'molecule_name': [],
                     'field': [],
                     'warning_string': []}

    # Check sdf filename & can be read
    validate_dict = check_sdf(sdf_file, validate_dict)

    suppl = Chem.SDMolSupplier(sdf_file)
    print('%d mols detected (including blank mol)' % (len(suppl),))
    blank_mol = suppl[0]
    if blank_mol is None:
        validate_dict = add_warning(molecule_name='Blank Mol',
                                    field='N/A',
                                    warning_string='your blank molecule could not be read by rdkit. The molecule must have at least one atom! No other checks were done',
                                    validate_dict=validate_dict)
        validated = False
        return validate_dict, validated
    validate_dict = check_compound_set(blank_mol, validate_dict)
    other_mols = []
    for i in range(1, len(suppl)):
        other_mols.append(suppl[i])

    # all mol checks
    # - all mols have the same properties
    all_props = []
    for mol in suppl:
        all_props.extend([key for key in mol.GetPropsAsDict().keys()])
    unique_props = list(set(all_props))
    for mol in suppl:
        props = [key for key in mol.GetPropsAsDict().keys()]
        diff_list = np.setdiff1d(props, unique_props)
        for diff in diff_list:
            add_warning(molecule_name=mol.GetProp('_Name'),
                        field='property (missing)',
                        warning_string='%s property is missing from this molecule' % (diff,),
                        validate_dict=validate_dict)

    # Check version in blank mol
    validate_dict = check_ver_name(blank_mol, version, validate_dict)

    # Check compuslory fields in blank mol props
    validate_dict = check_blank_mol_props(blank_mol, validate_dict)

    # Check properties have been described and validate url
    validate_dict = check_blank_prop(blank_mol, validate_dict)

    # main mols checks
    # - missing compulsary fields
    # - check name characters
    # - check pdb assignment and if pdb filepath exists
    # - check compulsory field populated
    # - check SMILES can be opended by rdkit
    # (check api for pdb if fragalysis)
    for m in other_mols:
        validate_dict = check_mol_props(m, validate_dict)
        validate_dict = check_name_characters(m.GetProp('_Name'), validate_dict)
        validate_dict = check_pdb(m, validate_dict, target, zfile)
        validate_dict = check_refmol(m, validate_dict, target)
        validate_dict = check_field_populated(m, validate_dict)
        validate_dict = check_SMILES(m, validate_dict)

    if len(validate_dict['molecule_name']) != 0:
        validated = False

    return validate_dict, validated
### End Validating ###