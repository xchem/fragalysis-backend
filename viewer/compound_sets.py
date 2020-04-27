import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
django.setup()

from rdkit import Chem
from viewer.models import (
    Target,
    CompoundSet,
    ComputedCompound,
    ScoreDescription,
    NumericalScoreValues,
    TextScoreValues,
    Protein,
    Molecule)
import ast
import os.path


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

def get_prot(mol, compound_set, filename):
    pdb_option = mol.GetProp('ref_pdb')
    path = os.path.abspath(filename)
    base = '/'.join(path.split('/')[:-1])
    if os.path.isdir(os.path.join(base, pdb_option)):
        field = os.path.join(base, pdb_option)
        print('PROT: ' + field)
    else:
        name = compound_set.target.title + '-' + pdb_option
        print('PROT: ' + name)
        prot = Protein.objects.get(code=name)
        field = prot.pdb_info

    return field


def set_props(cpd, props, compound_set):
    if 'ref_mols' and 'ref_pdb' not in props.keys():
        raise Exception('ref_mols and ref_pdb not set!')
    set_obj = ScoreDescription.objects.filter(compound_set=compound_set)
    set_props_list = [s.name for s in set_obj]
    for key in props.keys():
        if key in set_props_list not in ['ref_mols', 'ref_pdb']:
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


def set_mol(mol, compound_set, filename):
    smiles = Chem.MolToSmiles(mol)
    name = mol.GetProp('_Name')
    mol_block = Chem.MolToMolBlock(mol)

    insp = mol.GetProp('ref_mols')
    insp = insp.split(',')
    insp = [i.strip() for i in insp]
    insp_frags = [Molecule.objects.get(prot_id__code=str(compound_set.target.title + '-' + i)) for i in insp]

    prot_field = get_prot(mol, compound_set, filename)

    cpd = ComputedCompound()
    cpd.sdf_info = mol_block
    cpd.compound_set = compound_set
    cpd.name = name
    cpd.smiles = smiles
    cpd.pdb_info = prot_field

    cpd.save()

    [cpd.inspiration_frags.add(mol) for mol in insp_frags]

    cpd.save()

    return cpd


def process_mol(mol, compound_set, filename):
    cpd = set_mol(mol, compound_set, filename)
    other_props = mol.GetPropsAsDict()
    compound_set = set_props(cpd, other_props, compound_set)

    return compound_set


def set_descriptions(filename, compound_set):
    suppl = Chem.SDMolSupplier(str(filename))
    description_mol = suppl[0]
    description_dict = description_mol.GetPropsAsDict()
    version = description_mol.GetProp('_Name')
    compound_set.spec_version = version.split('_')[-1]
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


def process_compound_set(target, filename):
    print('processing compound set: ' + filename)
    filename = str(filename)
    # create a new compound set
    set_name = ''.join(filename.split('/')[-1].replace('.sdf','').split('_')[1:])
    compound_set = CompoundSet()
    compound_set.name = set_name
    matching_target = Target.objects.get(title=target)
    compound_set.target = matching_target
    compound_set.submitted_sdf = filename


    # set descriptions and get all other mols back
    mols_to_process = set_descriptions(filename=filename, compound_set=compound_set)

    # process every other mol
    for mol in mols_to_process:
        process_mol(mol, compound_set, filename)

    # check that molecules have been added to the compound set
    check = ComputedCompound.objects.filter(compound_set=compound_set)
    print(str(len(check)) + '/' + str(len(mols_to_process)) + ' succesfully processed in ' + set_name + ' cpd set')

    # save the compound set
    compound_set.save()

    # if no molecules were processed, delete the compound set
    if len(check) == 0:
        compound_set.delete()
        print('No molecules processed... deleting ' + set_name + ' compound set')

