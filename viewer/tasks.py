import os
import datetime
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
django.setup()
from django.conf import settings

from rdkit import Chem
from viewer.models import Target
import os.path

from celery import shared_task
from sdf_check import *
from compound_set_upload import *

# Bit to check if redis and celery services are working and available
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

### Uploading ###

@shared_task
def process_compound_set(validate_output):
    # Validate output is a tuple - this is one way to get
    # Celery chaining to work where second function uses tuple output
    # from first function called
    validate_dict, validated, target, filename, zfile = validate_output

    if not validated:
        return (validate_dict,validated)

    if validated:
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

@shared_task
def validate(sdf_file, target=None, zfile=None):
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

    return (validate_dict, validated, sdf_file, target, zfile)
### End Validating ###
