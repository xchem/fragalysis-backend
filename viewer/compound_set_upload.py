import os
import datetime
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
django.setup()
from django.conf import settings

from rdkit import Chem
from viewer.models import (
    ScoreDescription,
    Protein,
    Target,
    ComputedSetSubmitter)
import os.path

def get_inspiration_frags(cpd, compound_set):
    # Don't need...
    del cpd
    del compound_set


def process_pdb(pdb_code, target, zfile, zfile_hashvals):
    # pdb_fn = zfile[pdb_code].split('/')[-1]

    if zfile_hashvals:
        for key in zfile_hashvals.keys():
            if key == pdb_code:
                pdb_code = f'{pdb_code}-{zfile_hashvals[pdb_code]}'

    pdb_fp = zfile[pdb_code]
    pdb_fn = zfile[pdb_code].split('/')[-1]

    # Move and save the protein pdb from tmp to pdbs folder
    # pdb may have already been moved and Protein object created

    # prot_objs = Protein.objects.filter(code=pdb_code)
    # if len(prot_objs) > 0:
    #     raise Exception(f'Something went wrong with pdb zip upload: {[c.code for c in prot_objs]}')

    ## THIS BIT ISN'T MOVING THE FILES PROPERLY
    new_filename = settings.MEDIA_ROOT + 'pdbs/' + pdb_fn
    old_filename = settings.MEDIA_ROOT + pdb_fp
    os.renames(old_filename, new_filename)

    # Create Protein object
    prot = Protein()
    prot.code = pdb_code
    target_obj = Target.objects.get(title=target)
    prot.target_id = target_obj
    prot.pdb_info = 'pdbs/' + pdb_fn
    prot.save()

    # Get Protein object
    prot_obj = Protein.objects.get(code=pdb_code)

    return prot_obj

# use zfile object for pdb files uploaded in zip
def get_prot(mol, target, compound_set, zfile, zfile_hashvals=None):
    pdb_fn = mol.GetProp('ref_pdb').split('/')[-1]

    if zfile:
        pdb_code = pdb_fn.replace('.pdb','')
        prot_obj = process_pdb(pdb_code=pdb_code, target=target, zfile=zfile, zfile_hashvals=zfile_hashvals)
        field = prot_obj.pdb_info

    else:
        name = compound_set.target.title + '-' + pdb_fn
        prot_obj = Protein.objects.get(code__contains=name.split(':')[0].split('_')[0])
        field = prot_obj.pdb_info

    return field


def get_submission_info(description_mol):
    y_m_d = description_mol.GetProp('generation_date').split('-')

    submitter = ComputedSetSubmitter.objects.get_or_create(name=description_mol.GetProp('submitter_name'),
                                                           method=description_mol.GetProp('method'),
                                                           email=description_mol.GetProp('submitter_email'),
                                                           institution=description_mol.GetProp('submitter_institution'),
                                                           generation_date=datetime.date(int(y_m_d[0]), int(y_m_d[1]), int(y_m_d[2])))[0]

    return submitter


def get_additional_mols(filename, compound_set):
    suppl = Chem.SDMolSupplier(str(filename))
    mols = []

    for i in range(0, len(suppl)):
        mols.append(suppl[i])

    descriptions_list = list(
        set([item for sublist in [list(m.GetPropsAsDict().keys()) for m in mols] for item in sublist]))

    missing = []

    for desc in descriptions_list:
        existing = ScoreDescription.objects.filter(computed_set=compound_set, name=desc)
        if len(existing)==0 and desc not in ['original SMILES', 'ref_mols', 'ref_pdb']:
            missing.append(desc)

    if len(missing)>0:
        return f"Missing score descriptions for: {', '.join(missing)}, please re-upload"

    return mols
