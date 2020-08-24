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
    ComputedMolecule,
    ScoreDescription,
    NumericalScoreValues,
    TextScoreValues,
    Protein,
    Molecule,
    ComputedSetSubmitter)
import ast
import os.path

import os

import psutil

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
        if type(t) in [int, int, float, bool]:
            if t in [
                True, False, 'TRUE', 'FALSE', 'true', 'false', 'yes', 'no', 'YES', 'NO', 'Yes', 'No', "Y", "N", "y", "n"
            ]:
                return 'BIT'
            if type(t) is int or type(t) is int:
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
        if pdb_option in zfile:
            field = zfile[pdb_option]

    else:
        name = compound_set.target.title + '-' + pdb_option
        print(('PROT: ' + name))
        prot = Protein.objects.get(code__contains=name.split('_')[0])
        field = prot.pdb_info

    return field


def set_props(cpd, props, compound_set):
    if 'ref_mols' and 'ref_pdb' not in list(props.keys()):
        raise Exception('ref_mols and ref_pdb not set!')
    set_obj = ScoreDescription.objects.filter(computed_set=compound_set)
    set_props_list = [s.name for s in set_obj]
    for key in list(props.keys()):
        if key in set_props_list not in ['ref_mols', 'ref_pdb', 'original SMILES']:
            if dataType(str(props[key]))=='TEXT':
                score_value = TextScoreValues()
            else:
                score_value = NumericalScoreValues()
            score_value.score = ScoreDescription.objects.get(computed_set=compound_set,
                                                             name=key)
            score_value.value = props[key]
            score_value.compound = cpd
            score_value.save()

    return set_obj


def set_mol(mol, compound_set, filename, zfile=None):
    # zfile = {'zip_obj': zf, 'zf_list': zip_names}

    smiles = Chem.MolToSmiles(mol)
    inchi = Chem.inchi.MolToInchi(mol)
    from .tasks import create_mol
    name = mol.GetProp('_Name')
    long_inchi=None
    if len(inchi)>255:
        long_inchi = inchi
        inchi = inchi[0:254]

    ref_cpd = create_mol(inchi, name=name, long_inchi=long_inchi)

    mol_block = Chem.MolToMolBlock(mol)

    insp = mol.GetProp('ref_mols')
    insp = insp.split(',')
    insp = [i.strip() for i in insp]
    insp_frags = []
    for i in insp:
        mols = Molecule.objects.filter(prot_id__code__contains=str(compound_set.target.title + '-' + i.split('_')[0]),
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
        # Save the compound set
        filename = settings.MEDIA_ROOT + prot_field
        prot_field = filename
        compound_set.save()

    #  need to add Compound before saving
    # see if anything exists already
    existing = ComputedMolecule.objects.filter(name=name, smiles=smiles, computed_set=compound_set)

    if len(existing)==1:
        cpd = existing[0]
    if len(existing)>1:
        [c.delete for c in existing]
        cpd = ComputedMolecule()
    else:
        cpd = ComputedMolecule()
    cpd.compound = ref_cpd
    cpd.computed_set = compound_set
    cpd.sdf_info = mol_block
    cpd.name = name
    cpd.smiles = smiles
    cpd.pdb_info = prot_field
    #cpd.original_smiles = orig
    cpd.save()

    [cpd.computed_inspirations.add(mol) for mol in insp_frags]

    cpd.save()

    return cpd


def process_mol(mol, compound_set, filename, zfile=None):
    cpd = set_mol(mol, compound_set, filename, zfile)
    other_props = mol.GetPropsAsDict()
    compound_set = set_props(cpd, other_props, compound_set)

    return compound_set


def get_submission_info(description_mol):
    y_m_d = description_mol.GetProp('generation_date').split('-')

    submitter = ComputedSetSubmitter.objects.get_or_create(name=description_mol.GetProp('submitter_name'),
                                                           method=description_mol.GetProp('method'),
                                                           email=description_mol.GetProp('submitter_email'),
                                                           institution=description_mol.GetProp('submitter_institution')
                                                           )[0]

    submitter.generation_date = datetime.date(int(y_m_d[0]), int(y_m_d[1]), int(y_m_d[2]))

    submitter.save()

    return submitter


def set_descriptions(filename, compound_set):
    # remove any old scores
    text_scores = TextScoreValues.objects.filter(score__computed_set=compound_set)
    [t.delete() for t in text_scores]
    num_scores = NumericalScoreValues.objects.filter(score__computed_set=compound_set)
    [n.delete() for n in num_scores]

    suppl = Chem.SDMolSupplier(str(filename))
    description_mol = suppl[0]

    mols = []

    for i in range(1, len(suppl)):
        mols.append(suppl[i])
    #
    # description_keys = []
    #
    # for i in range(1, len(suppl)):
    #     description_keys.extend(suppl[i].GetPropsAsDict().keys())

    descriptions_needed = list(set([item for sublist in [list(m.GetPropsAsDict().keys()) for m in mols] for item in sublist]))
     # list(set())

    submitter = get_submission_info(description_mol)

    description_dict = description_mol.GetPropsAsDict()
    version = description_mol.GetProp('_Name')
    compound_set.spec_version = version.split('_')[-1]
    method = description_mol.GetProp('ref_url')
    compound_set.method_url = method
    compound_set.submitter = submitter
    compound_set.save()

    for key in list(description_dict.keys()):
        if key in descriptions_needed and key not in ['ref_mols', 'ref_pdb', 'index', 'Name', 'original SMILES']:
            desc = ScoreDescription.objects.get_or_create(computed_set=compound_set,
                                                          name=key,
                                                          description=description_dict[key],
                                                          )[0]

    return mols

