import ast
import datetime
import logging
import os
import shutil
import uuid
import zipfile
from typing import Any, Dict, List, Optional, Tuple

from openpyxl.utils import get_column_letter

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django

django.setup()
from django.conf import settings

logger = logging.getLogger(__name__)

from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from viewer.models import (
    Compound,
    ComputedMolecule,
    ComputedSet,
    ComputedSetSubmitter,
    NumericalScoreValues,
    ScoreDescription,
    SiteObservation,
    Target,
    TextScoreValues,
    User,
)
from viewer.utils import add_props_to_sdf_molecule, is_url, word_count


def dataType(a_str: str) -> str:
    lean_str = a_str.strip()
    if not lean_str:
        return 'BLANK'

    try:
        t = ast.literal_eval(lean_str)
    except (ValueError, SyntaxError):
        return 'TEXT'
    else:
        if type(t) in [int, int, float, bool]:
            if t in [
                True,
                False,
                'TRUE',
                'FALSE',
                'true',
                'false',
                'yes',
                'no',
                'YES',
                'NO',
                'Yes',
                'No',
                "Y",
                "N",
                "y",
                "n",
            ]:
                return 'BIT'
            if type(t) is int or type(t) is int:
                return 'INT'
            if type(t) is float:
                return 'FLOAT'

            # Can't get here?
            assert False
        else:
            return 'TEXT'


class PdbOps:
    def save_pdb_zip(
        self, pdb_file
    ) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, str]]]:
        zfile = None
        zfile_hashvals = None
        if pdb_file:
            zf = zipfile.ZipFile(pdb_file)
            zip_lst = zf.namelist()
            zfile = {}
            zfile_hashvals = {}
            for filename in zip_lst:
                # only handle pdb files
                if filename.split('.')[-1] == 'pdb':
                    # Test if Protein object already exists
                    code = filename.split('/')[-1].replace('.pdb', '')
                    test_pdb_code = filename.split('/')[-1].replace('.pdb', '')
                    test_site_obvs_objs = SiteObservation.objects.filter(
                        code=test_pdb_code
                    )

                    if len(test_site_obvs_objs) != 0:
                        # make a unique pdb code as not to overwrite existing object
                        rand_str = uuid.uuid4().hex
                        test_pdb_code = f'{code}#{rand_str}'
                        zfile_hashvals[code] = rand_str

                    fn = f'{test_pdb_code}.pdb'
                    pdb_path = default_storage.save(
                        f'tmp/{fn}', ContentFile(zf.read(filename))
                    )
                    zfile[test_pdb_code] = pdb_path

            # Close the zip file
            if zf:
                zf.close()

        return zfile, zfile_hashvals

    def run(self, params) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, str]]]:
        return self.save_pdb_zip(params['pdb_zip'])


class MolOps:
    def __init__(
        self,
        user_id,
        sdf_filename,
        submitter_name,
        submitter_method,
        target,
        version,
        zfile,
        zfile_hashvals,
    ):
        self.user_id = user_id
        self.sdf_filename = sdf_filename
        self.submitter_name = submitter_name
        self.submitter_method = submitter_method
        self.target = target
        self.version = version
        self.zfile = zfile
        self.zfile_hashvals = zfile_hashvals

    def process_pdb(self, pdb_code, target, zfile, zfile_hashvals) -> SiteObservation:
        for key in zfile_hashvals.keys():
            if key == pdb_code:
                pdb_code = f'{pdb_code}#{zfile_hashvals[pdb_code]}'

        pdb_fp = zfile[pdb_code]
        pdb_fn = zfile[pdb_code].split('/')[-1]

        new_filename = f'{settings.MEDIA_ROOT}pdbs/{pdb_fn}'
        old_filename = settings.MEDIA_ROOT + pdb_fp
        shutil.copy(old_filename, new_filename)

        # Create Protein object
        target_obj = Target.objects.get(title=target)
        # prot.target_id = target_obj
        site_obvs, created = SiteObservation.objects.get_or_create(
            code=pdb_code, target_id=target_obj
        )
        # prot.code = pdb_code
        if created:
            target_obj = Target.objects.get(title=target)
            site_obvs.target_id = target_obj
            site_obvs.pdb_info = f'pdbs/{pdb_fn}'
            site_obvs.save()

        return site_obvs

    # use zfile object for pdb files uploaded in zip
    def get_site_observation(
        self, property_name, mol, target, compound_set, zfile, zfile_hashvals
    ) -> Optional[SiteObservation]:
        # Get a SiteObservation from the molecule using
        # a named property (i.e. lhs_pdb or ref_pdb for example)

        if not mol.HasProp(property_name):
            logger.warning(
                'Molecule %s has no "%s" property (%s, %s)',
                mol,
                property_name,
                target,
                compound_set,
            )
            return None

        pdb_fn = mol.GetProp(property_name).split('/')[-1]
        site_obvs = None

        if zfile:
            pdb_code = pdb_fn.replace('.pdb', '')
            site_obvs = self.process_pdb(
                pdb_code=pdb_code,
                target=target,
                zfile=zfile,
                zfile_hashvals=zfile_hashvals,
            )
        else:
            name = pdb_fn
            try:
                site_obvs = SiteObservation.objects.get(
                    code__contains=name,
                    experiment__experiment_upload__target__title=target,
                )
            except SiteObservation.DoesNotExist:
                # Initial SiteObservation lookup failed.
                logger.warning(
                    'Failed to get SiteObservation object (target=%s name=%s)',
                    compound_set.target.title,
                    name,
                )
                # Try alternatives.
                # If all else fails then the site_obvs will be 'None'
                qs = SiteObservation.objects.filter(
                    code__contains=name,
                    experiment__experiment_upload__target__title=target,
                )
                if qs.exists():
                    logger.info(
                        'Found SiteObservation containing name=%s qs=%s',
                        name,
                        qs,
                    )
                else:
                    alt_name = name.split(':')[0].split('_')[0]
                    qs = SiteObservation.objects.filter(
                        code__contains=alt_name,
                        experiment__experiment_upload__target__title=target,
                    )
                    if qs.exists():
                        logger.info(
                            'Found SiteObservation containing alternative name=%s qs=%s',
                            alt_name,
                            qs,
                        )
                if qs.count() > 0:
                    logger.debug(
                        'Found alternative (target=%s name=%s)',
                        compound_set.target.title,
                        name,
                    )
                    site_obvs = qs[0]

        if not site_obvs:
            logger.warning(
                'No SiteObservation found (target=%s pdb_fn=%s)',
                compound_set.target.title,
                pdb_fn,
            )

        return site_obvs

    def create_mol(self, inchi, long_inchi=None, name=None) -> Compound:
        # check for an existing compound, returning a Compound
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

    def set_props(self, cpd, props, compound_set) -> List[ScoreDescription]:
        if 'ref_mols' and 'ref_pdb' not in list(props.keys()):
            raise Exception('ref_mols and ref_pdb not set!')
        set_obj = ScoreDescription.objects.filter(computed_set=compound_set)
        assert set_obj

        set_props_list = [s.name for s in set_obj]
        for key in list(props.keys()):
            if key in set_props_list not in ['ref_mols', 'ref_pdb', 'original SMILES']:
                if dataType(str(props[key])) == 'TEXT':
                    score_value = TextScoreValues()
                else:
                    score_value = NumericalScoreValues()
                score_value.score = ScoreDescription.objects.get(
                    computed_set=compound_set, name=key
                )
                score_value.value = props[key]
                score_value.compound = cpd
                score_value.save()

        return set_obj

    def set_mol(
        self, mol, target, compound_set, filename, zfile=None, zfile_hashvals=None
    ) -> ComputedMolecule:
        # Don't need...
        assert target
        assert compound_set

        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.inchi.MolToInchi(mol)
        molecule_name = mol.GetProp('_Name')
        long_inchi = None
        if len(inchi) > 255:
            long_inchi = inchi
            inchi = inchi[:254]

        compound: Compound = self.create_mol(
            inchi, name=molecule_name, long_inchi=long_inchi
        )

        insp = mol.GetProp('ref_mols')
        insp = insp.split(',')
        insp = [i.strip() for i in insp]
        insp_frags = []
        for i in insp:
            # try exact match first
            try:
                site_obvs = SiteObservation.objects.get(
                    code=str(i),
                    experiment__experiment_upload__target_id=compound_set.target,
                )
                ref = site_obvs
            except SiteObservation.DoesNotExist:
                qs = SiteObservation.objects.filter(
                    code=str(i.split(':')[0].split('_')[0]),
                    experiment__experiment_upload__target_id=compound_set.target,
                )
                if not qs.exists():
                    raise Exception(  # pylint: disable=raise-missing-from
                        'No matching molecules found for inspiration frag ' + i
                    )

                if qs.count() > 1:
                    ids = [m.cmpd.id for m in qs]
                    ind = ids.index(max(ids))
                    ref = qs[ind]
                elif qs.count() == 1:
                    ref = qs[0]

            insp_frags.append(ref)

        # Try to get the LHS SiteObservation,
        # This will be used to set the ComputedMolecule.site_observation_code.
        # This may fail.
        lhs_property = 'lhs_pdb'
        lhs_so = self.get_site_observation(
            lhs_property,
            mol,
            target,
            compound_set,
            zfile,
            zfile_hashvals=zfile_hashvals,
        )
        if not lhs_so:
            logger.warning(
                'Failed to get a LHS SiteObservation (%s) for %s, %s, %s',
                lhs_property,
                mol,
                target,
                compound_set,
            )

        # Try to get the reference SiteObservation,
        # This will be used to set the ComputedMolecule.reference_code.
        # This may fail.
        ref_property = 'ref_pdb'
        ref_so = self.get_site_observation(
            ref_property,
            mol,
            target,
            compound_set,
            zfile,
            zfile_hashvals=zfile_hashvals,
        )
        if not ref_so:
            logger.warning(
                'Failed to get a Reference SiteObservation (%s) for %s, %s, %s',
                ref_property,
                mol,
                target,
                compound_set,
            )

        # A LHS or Reference protein must be provided.
        # (Part of "Fix behaviour of RHS [P] button - also RHS upload change", issue #1249)
        if not lhs_so and not ref_so:
            logger.error(
                'ComputedMolecule has no LHS (%s) or Reference (%s) property',
                lhs_property,
                ref_property,
            )

        # Need a ComputedMolecule before saving.
        # Check if anything exists already...
        existing_computed_molecules = ComputedMolecule.objects.filter(
            molecule_name=molecule_name, smiles=smiles, computed_set=compound_set
        )

        computed_molecule: Optional[ComputedMolecule] = None
        if len(existing_computed_molecules) == 1:
            logger.info(
                'Using existing ComputedMolecule %s', existing_computed_molecules[0]
            )
            computed_molecule = existing_computed_molecules[0]
        elif len(existing_computed_molecules) > 1:
            logger.warning('Deleting existing ComputedMolecules (more than 1 found')
            for exist in existing_computed_molecules:
                logger.info('Deleting ComputedMolecule %s', exist)
                exist.delete()
            computed_molecule = ComputedMolecule()
        if not computed_molecule:
            logger.info('Creating new ComputedMolecule')
            computed_molecule = ComputedMolecule()

        assert computed_molecule
        computed_molecule.compound = compound
        computed_molecule.computed_set = compound_set
        computed_molecule.sdf_info = Chem.MolToMolBlock(mol)
        computed_molecule.site_observation_code = lhs_so.code if lhs_so else None
        computed_molecule.reference_code = ref_so.code if ref_so else None
        computed_molecule.molecule_name = molecule_name
        computed_molecule.name = f"{target}-{computed_molecule.identifier}"
        computed_molecule.smiles = smiles
        # Extract possible reference URL and Rationale
        # URLs have to be valid URLs and rationals must contain more than one word
        ref_url: Optional[str] = (
            mol.GetProp('ref_url') if mol.HasProp('ref_url') else None
        )
        computed_molecule.ref_url = ref_url if is_url(ref_url) else None
        rationale: Optional[str] = (
            mol.GetProp('rationale') if mol.HasProp('rationale') else None
        )
        computed_molecule.rationale = rationale if word_count(rationale) > 1 else None
        # To avoid the error...
        #   needs to have a value for field "id"
        #   before this many-to-many relationship can be used.
        # We must save this ComputedMolecule to generate an "id"
        # before adding inspirations
        computed_molecule.save()
        for insp_frag in insp_frags:
            computed_molecule.computed_inspirations.add(insp_frag)
        # Done
        computed_molecule.save()

        # No update the molecule in the original file...
        add_props_to_sdf_molecule(
            sdf_file=filename,
            molecule=molecule_name,
            properties={"target_identifier": computed_molecule.name},
        )

        return computed_molecule

    def get_submission_info(self, description_mol) -> ComputedSetSubmitter:
        y_m_d = description_mol.GetProp('generation_date').split('-')
        return ComputedSetSubmitter.objects.get_or_create(
            name=description_mol.GetProp('submitter_name'),
            method=description_mol.GetProp('method'),
            email=description_mol.GetProp('submitter_email'),
            institution=description_mol.GetProp('submitter_institution'),
            generation_date=datetime.date(int(y_m_d[0]), int(y_m_d[1]), int(y_m_d[2])),
        )[0]

    def process_mol(
        self, mol, target, compound_set, filename, zfile=None, zfile_hashvals=None
    ) -> List[ScoreDescription]:
        cpd = self.set_mol(mol, target, compound_set, filename, zfile, zfile_hashvals)
        other_props = mol.GetPropsAsDict()
        return self.set_props(cpd, other_props, compound_set)

    def set_descriptions(
        self, filename, computed_set: ComputedSet
    ) -> List[Chem.rdchem.Mol]:
        suppl = Chem.SDMolSupplier(str(filename))
        description_mol = suppl[0]

        mols = [suppl[i] for i in range(1, len(suppl))]
        descriptions_needed = list(
            {
                item
                for sublist in [list(m.GetPropsAsDict().keys()) for m in mols]
                for item in sublist
            }
        )

        computed_set.submitter = self.get_submission_info(description_mol)
        if description_mol.HasProp('ref_url'):
            computed_set.method_url = description_mol.GetProp('ref_url')
        computed_set.save()

        description_dict = description_mol.GetPropsAsDict()
        for key in description_dict.keys():
            if key in descriptions_needed and key not in [
                'ref_mols',
                'ref_pdb',
                'index',
                'Name',
                'original SMILES',
            ]:
                _ = ScoreDescription.objects.get_or_create(
                    computed_set=computed_set,
                    name=key,
                    description=description_dict[key],
                )

        return mols

    def task(self) -> ComputedSet:
        # Truncate submitted method (lower-case)?
        truncated_submitter_method: str = 'unspecified'
        if self.submitter_method:
            truncated_submitter_method = self.submitter_method[
                : ComputedSet.LENGTH_METHOD_IN_NAME
            ]
            if len(self.submitter_method) > len(truncated_submitter_method):
                logger.warning(
                    'ComputedSet submitter method is too long (%s). Truncated to "%s"',
                    self.submitter_method,
                    truncated_submitter_method,
                )
        else:
            logger.warning(
                'ComputedSet submitter method is not set. Using "%s"',
                truncated_submitter_method,
            )

        # Do we have any existing ComputedSets?
        # Ones with the same method and upload date?
        today: datetime.date = datetime.date.today()
        existing_sets: List[ComputedSet] = ComputedSet.objects.filter(
            method=truncated_submitter_method, upload_date=today
        ).all()
        # If so, find the one with the highest ordinal.
        latest_ordinal: int = 0
        for exiting_set in existing_sets:
            assert exiting_set.md_ordinal > 0
            if exiting_set.md_ordinal > latest_ordinal:
                latest_ordinal = exiting_set.md_ordinal
        if latest_ordinal:
            logger.info(
                'Found existing ComputedSets for method "%s" on %s (%d) with ordinal=%d',
                truncated_submitter_method,
                str(today),
                len(existing_sets),
                latest_ordinal,
            )
        # ordinals are 1-based
        new_ordinal: int = latest_ordinal + 1

        # The computed set "name" consists of the "method",
        # today's date and a 2-digit ordinal. The ordinal
        # is used to distinguish between computed sets uploaded
        # with the same method on the same day.
        assert new_ordinal > 0
        cs_name: str = f"{truncated_submitter_method}-{str(today)}-{get_column_letter(new_ordinal)}"
        logger.info('Creating new ComputedSet "%s"', cs_name)

        computed_set: ComputedSet = ComputedSet()
        computed_set.name = cs_name
        computed_set.md_ordinal = new_ordinal
        computed_set.upload_date = today
        computed_set.method = self.submitter_method[: ComputedSet.LENGTH_METHOD]
        computed_set.target = Target.objects.get(title=self.target)
        computed_set.spec_version = float(self.version.strip('ver_'))
        if self.user_id:
            computed_set.owner_user = User.objects.get(id=self.user_id)
        else:
            # The User ID may only be None if AUTHENTICATE_UPLOAD is False.
            # Here the ComputedSet owner will take on a default (anonymous) value.
            assert settings.AUTHENTICATE_UPLOAD is False
        computed_set.save()

        # Set descriptions in return for the Molecules.
        # This also sets the submitter and method URL properties of the computed set
        # while also saving it.
        sdf_filename = str(self.sdf_filename)
        mols_to_process = self.set_descriptions(
            filename=sdf_filename, computed_set=computed_set
        )

        # Process the molecules
        logger.info('%s mols_to_process=%s', computed_set, len(mols_to_process))
        for i in range(len(mols_to_process)):
            _ = self.process_mol(
                mols_to_process[i],
                self.target,
                computed_set,
                sdf_filename,
                self.zfile,
                self.zfile_hashvals,
            )

        # check compound set folder exists.
        cmp_set_folder = os.path.join(
            settings.MEDIA_ROOT, settings.COMPUTED_SET_MEDIA_DIRECTORY
        )
        if not os.path.isdir(cmp_set_folder):
            logger.info('Making ComputedSet folder (%s)', cmp_set_folder)
            os.mkdir(cmp_set_folder)

        # move and save the compound set
        new_filename = f'{settings.MEDIA_ROOT}{settings.COMPUTED_SET_MEDIA_DIRECTORY}/{computed_set.name}.sdf'
        os.rename(sdf_filename, new_filename)
        computed_set.submitted_sdf = sdf_filename
        computed_set.written_sdf_filename = new_filename
        computed_set.save()

        logger.info('Created %s', computed_set)

        return computed_set


def blank_mol_vals(sdf_file) -> Tuple[str, str, str]:
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
