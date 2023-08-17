import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
django.setup()
from django.conf import settings

import zipfile
import datetime
import ast
import shutil


# import the logging library
import logging
# Get an instance of a logger
logger = logging.getLogger(__name__)

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
import uuid

from viewer.models import (
    Compound,
    ComputedMolecule,
    ComputedSet,
    ComputedSetSubmitter,
    NumericalScoreValues,
    SiteObservation,
    ScoreDescription,
    Target,
    TextScoreValues,
    User,
)


def dataType(a_str):
    lean_str = a_str.strip()
    if len(lean_str) == 0: return 'BLANK'
    try:
        t = ast.literal_eval(lean_str)

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


class PdbOps:

    def save_pdb_zip(self, pdb_file):
        zfile = None
        zfile_hashvals = None
        if pdb_file:
            zf = zipfile.ZipFile(pdb_file)
            zip_lst = zf.namelist()
            zfile = {}
            zfile_hashvals = {}
            print(zip_lst)
            for filename in zip_lst:
                # only handle pdb files
                if filename.split('.')[-1] == 'pdb':
                    # Test if Protein object already exists
                    code = filename.split('/')[-1].replace('.pdb', '')
                    test_pdb_code = filename.split('/')[-1].replace('.pdb', '')
                    test_site_obvs_objs = SiteObservation.objects.filter(code=test_pdb_code)
                    print([c.code for c in test_site_obvs_objs])

                    if len(test_site_obvs_objs) != 0:
                        # make a unique pdb code as not to overwrite existing object
                        rand_str = uuid.uuid4().hex
                        test_pdb_code = f'{code}#{rand_str}'
                        zfile_hashvals[code] = rand_str

                    fn = test_pdb_code + '.pdb'

                    pdb_path = default_storage.save('tmp/' + fn,
                                                    ContentFile(zf.read(filename)))
                    zfile[test_pdb_code] = pdb_path

            # Close the zip file
            if zf:
                zf.close()

        return zfile, zfile_hashvals


    def run(self, params):
        zfile, zfile_hashval = self.save_pdb_zip(params['pdb_zip'])

        return zfile, zfile_hashval


class MolOps:

    def __init__(self, user_id, sdf_filename, submitter_name, submitter_method, target, version, zfile, zfile_hashvals):
        self.user_id = user_id
        self.sdf_filename = sdf_filename
        self.submitter_name = submitter_name
        self.submitter_method = submitter_method
        self.target = target
        self.version = version
        self.zfile = zfile
        self.zfile_hashvals = zfile_hashvals

    def process_pdb(self, pdb_code, target, zfile, zfile_hashvals):

        for key in zfile_hashvals.keys():
            if key == pdb_code:
                pdb_code = f'{pdb_code}#{zfile_hashvals[pdb_code]}'

        pdb_fp = zfile[pdb_code]
        pdb_fn = zfile[pdb_code].split('/')[-1]

        new_filename = settings.MEDIA_ROOT + 'pdbs/' + pdb_fn
        old_filename = settings.MEDIA_ROOT + pdb_fp
        shutil.copy(old_filename, new_filename)

        # Create Protein object
        target_obj = Target.objects.get(title=target)
        # prot.target_id = target_obj
        site_obvs, created = SiteObservation.objects.get_or_create(code=pdb_code, target_id=target_obj)
        # prot.code = pdb_code
        if created:
            target_obj = Target.objects.get(title=target)
            site_obvs.target_id = target_obj
            site_obvs.pdb_info = 'pdbs/' + pdb_fn
            site_obvs.save()

        return site_obvs

    # use zfile object for pdb files uploaded in zip
    def get_prot(self, mol, target, compound_set, zfile, zfile_hashvals):
        # The returned protein object may be None

        pdb_fn = mol.GetProp('ref_pdb').split('/')[-1]
        site_obvs = None

        if zfile:
            pdb_code = pdb_fn.replace('.pdb', '')
            site_obvs = self.process_pdb(pdb_code=pdb_code, target=target, zfile=zfile, zfile_hashvals=zfile_hashvals)
        else:
            name = compound_set.target.title + '-' + pdb_fn

            # try to get single exact match
            # name.split(':')[0].split('_')[0]
            try:
                site_obvs = SiteObservation.objects.get(code__contains=name)
            except SiteObservation.DoesNotExist:
                # Protein lookup failed.
                logger.warning('Failed to get Protein object (target=%s name=%s)',
                               compound_set.target.title, name)
                # Try an alternative.
                # If all else fails then the prot_obj will be 'None'
                qs = SiteObservation.objects.filter(code__contains=name)
                if not qs.exists():
                    qs = SiteObservation.objects.filter(code__contains=name.split(':')[0].split('_')[0])
                if qs.count() > 0:
                    site_obvs = qs[0]

        if not site_obvs:
            logger.warning('No Protein object (target=%s pdb_fn=%s)',
                           compound_set.target.title, pdb_fn)

        return site_obvs

    def create_mol(self, inchi, long_inchi=None, name=None):
        # check for an existing compound
        if long_inchi:
            cpd = Compound.objects.filter(long_inchi=long_inchi)
            sanitized_mol = Chem.MolFromInchi(long_inchi, sanitize=True)
        else:
            cpd = Compound.objects.filter(inchi=inchi)
            sanitized_mol = Chem.MolFromInchi(inchi, sanitize=True)

        if len(cpd) != 0:
            new_mol = cpd[0]
        elif len(cpd) == 0:
            # add molecule and return the object
            new_mol = Compound()

        new_mol.smiles = Chem.MolToSmiles(sanitized_mol)
        new_mol.inchi = inchi
        if long_inchi:
            new_mol.long_inchi = long_inchi
        new_mol.identifier = name

        # descriptors
        new_mol.mol_log_p = Crippen.MolLogP(sanitized_mol)
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


    def set_props(self, cpd, props, compound_set):
        if 'ref_mols' and 'ref_pdb' not in list(props.keys()):
            raise Exception('ref_mols and ref_pdb not set!')
        set_obj = ScoreDescription.objects.filter(computed_set=compound_set)
        set_props_list = [s.name for s in set_obj]
        for key in list(props.keys()):
            if key in set_props_list not in ['ref_mols', 'ref_pdb', 'original SMILES']:
                if dataType(str(props[key])) == 'TEXT':
                    score_value = TextScoreValues()
                else:
                    score_value = NumericalScoreValues()
                score_value.score = ScoreDescription.objects.get(computed_set=compound_set,
                                                                 name=key)
                score_value.value = props[key]
                score_value.compound = cpd
                score_value.save()

        return set_obj


    def set_mol(self, mol, target, compound_set, filename, zfile=None, zfile_hashvals=None):
        # Don't need...
        del filename

        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.inchi.MolToInchi(mol)
        name = mol.GetProp('_Name')
        long_inchi = None
        if len(inchi) > 255:
            long_inchi = inchi
            inchi = inchi[:254]

        ref_cpd = self.create_mol(inchi, name=name, long_inchi=long_inchi)

        mol_block = Chem.MolToMolBlock(mol)

        insp = mol.GetProp('ref_mols')
        insp = insp.split(',')
        insp = [i.strip() for i in insp]
        insp_frags = []
        for i in insp:
            # try exact match first
            try:
                site_obvs = SiteObservation.objects.get(
                    prot_id__code__contains=str(compound_set.target.title + '-' + i),
                    prot_id__target_id=compound_set.target,
                )
                ref = site_obvs
            except SiteObservation.DoesNotExist:

                qs = SiteObservation.objects.filter(
                    prot_id__code__contains=str(compound_set.target.title + '-' + i.split(':')[0].split('_')[0]),
                    prot_id__target_id=compound_set.target,
                )
                if not qs.exists():
                    raise Exception('No matching molecules found for inspiration frag ' + i)  # pylint: disable=raise-missing-from

                if qs.count() > 1:
                    ids = [m.cmpd_id.id for m in qs]
                    ind = ids.index(max(ids))
                    ref = qs[ind]
                elif qs.count()  == 1:
                    ref = qs[0]

            insp_frags.append(ref)

        _ = mol.GetProp('original SMILES')

        # Try to get the protein object.
        # This may fail.
        prot = self.get_prot(mol, target, compound_set, zfile, zfile_hashvals=zfile_hashvals)
        if not prot:
            logger.warning('get_prot() failed to return a Protein object')

        #  need to add Compound before saving
        # see if anything exists already
        existing = ComputedMolecule.objects.filter(name=name, smiles=smiles, computed_set=compound_set)

        if len(existing) == 1:
            cpd = existing[0]
        if len(existing) > 1:
            for exist in existing:
                exist.delete()
            cpd = ComputedMolecule()
        elif len(existing) == 0:
            cpd = ComputedMolecule()

        cpd.compound = ref_cpd
        cpd.computed_set = compound_set
        cpd.sdf_info = mol_block
        cpd.name = name
        cpd.smiles = smiles
        cpd.pdb = prot
        cpd.save()

        for insp_frag in insp_frags:
            cpd.computed_inspirations.add(insp_frag)

        cpd.save()

        return cpd


    def get_submission_info(self, description_mol):
        y_m_d = description_mol.GetProp('generation_date').split('-')

        submitter = ComputedSetSubmitter.objects.get_or_create(name=description_mol.GetProp('submitter_name'),
                                                               method=description_mol.GetProp('method'),
                                                               email=description_mol.GetProp('submitter_email'),
                                                               institution=description_mol.GetProp(
                                                                   'submitter_institution'),
                                                               generation_date=datetime.date(int(y_m_d[0]),
                                                                                             int(y_m_d[1]),
                                                                                             int(y_m_d[2])))[0]

        return submitter

    def process_mol(self, mol, target, compound_set, filename, zfile=None, zfile_hashvals=None):
        cpd = self.set_mol(mol, target, compound_set, filename, zfile, zfile_hashvals)
        other_props = mol.GetPropsAsDict()
        compound_set = self.set_props(cpd, other_props, compound_set)

        return compound_set

    def set_descriptions(self, filename, compound_set):

        suppl = Chem.SDMolSupplier(str(filename))
        description_mol = suppl[0]

        mols = []

        for i in range(1, len(suppl)):
            mols.append(suppl[i])

        descriptions_needed = list(
            set([item for sublist in [list(m.GetPropsAsDict().keys()) for m in mols] for item in sublist]))

        submitter = self.get_submission_info(description_mol)

        description_dict = description_mol.GetPropsAsDict()
        version = description_mol.GetProp('_Name')
        compound_set.spec_version = version.split('_')[-1]
        method = description_mol.GetProp('ref_url')
        compound_set.method_url = method
        compound_set.submitter = submitter
        compound_set.save()

        for key in list(description_dict.keys()):
            if key in descriptions_needed and key not in ['ref_mols', 'ref_pdb', 'index', 'Name', 'original SMILES']:
                _ = ScoreDescription.objects.get_or_create(computed_set=compound_set,
                                                           name=key,
                                                           description=description_dict[key],
                                                           )

        return mols

    def task(self):
        user = User.objects.get(id=self.user_id)
        sdf_filename = str(self.sdf_filename)

        # create a new compound set
        set_name = ''.join(sdf_filename.split('/')[-1].replace('.sdf', '').split('_')[1:])
        unique_name = "".join(self.submitter_name.split()) + '-' + "".join(self.submitter_method.split())

        existing = ComputedSet.objects.filter(unique_name=unique_name)
        len_existing = len(existing)
        logger.info('Existing ComputedSet.unique_name=%s len(existing)=%s',
                    unique_name, len_existing)

        if len_existing == 1:
            logger.info('Using existing ComputedSet')
            compound_set = existing[0]
        elif len_existing > 1:
            raise Exception('Too many ComputedSet instances exist'
                            f' (unique_name="{unique_name}" len_existing={len_existing})')
        else:
            compound_set = ComputedSet()

        text_scores = TextScoreValues.objects.filter(score__computed_set=compound_set)
        num_scores = NumericalScoreValues.objects.filter(score__computed_set=compound_set)

        old_mols = [o.compound for o in text_scores]
        old_mols.extend([o.compound for o in num_scores])
        print(list(set(old_mols)))

        compound_set.name = set_name
        matching_target = Target.objects.get(title=self.target)
        compound_set.target = matching_target
        ver = float(self.version.strip('ver_'))
        compound_set.spec_version = ver
        compound_set.unique_name = "".join(self.submitter_name.split()) + '-' + "".join(self.submitter_method.split())
        compound_set.owner_user = user
        compound_set.save()

        # set descriptions and get all other mols back
        mols_to_process = self.set_descriptions(filename=sdf_filename, compound_set=compound_set)

        # process every other mol
        for i in range(0, len(mols_to_process)):
            self.process_mol(mols_to_process[i], self.target, compound_set, sdf_filename, self.zfile, self.zfile_hashvals)

        # check that molecules have been added to the compound set
        _ = ComputedMolecule.objects.filter(computed_set=compound_set)

        # check compound set folder exists.
        cmp_set_folder = os.path.join(settings.MEDIA_ROOT, 'compound_sets')
        if not os.path.isdir(cmp_set_folder):
            os.mkdir(cmp_set_folder)

        # move and save the compound set
        new_filename = settings.MEDIA_ROOT + 'compound_sets/' + sdf_filename.split('/')[-1]
        shutil.copy(sdf_filename, new_filename)
        # os.renames(sdf_filename, new_filename)
        compound_set.submitted_sdf = new_filename
        compound_set.save()

        # old_mols = [o.compound for o in old_s2]
        # old_mols.extend([o.compound for o in old_s1])
        # print(list(set(old_mols)))

        for old_mol in old_mols:
            old_mol.delete()

        return compound_set


def blank_mol_vals(sdf_file):
    """Returns the submitter name, method and version (_Name) if present.
    If not present the corresponding values are empty strings.
    """
    suppl = Chem.SDMolSupplier(sdf_file)
    if not suppl:
        return '', '', ''
    # print('%d mols detected (including blank mol)' % (len(suppl),))
    blank_mol = suppl[0]
    if not blank_mol:
        return '', '', ''

    # Get submitter name/info for passing into upload to get unique name
    submitter_name = ''
    if blank_mol.HasProp('submitter_name'):
        submitter_name = blank_mol.GetProp('submitter_name')

    submitter_method = ''
    if blank_mol.HasProp('method'):
        submitter_method = blank_mol.GetProp('method')

    version = ''
    if blank_mol.HasProp('_Name'):
        version = blank_mol.GetProp('_Name')

    return submitter_name, submitter_method, version
