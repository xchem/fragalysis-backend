import ast
import datetime
import logging
import os
import uuid
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from openpyxl.utils import get_column_letter

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django

django.setup()

from django.conf import settings
from django.core.exceptions import MultipleObjectsReturned
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.db.models import F, TextField, Value
from django.db.models.expressions import Func
from rdkit import Chem

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
from viewer.utils import add_props_to_sdf_molecule, alphanumerator, is_url, word_count

logger = logging.getLogger(__name__)


# maximum distance between corresponding atoms in poses
_DIST_LIMIT = 0.5


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
        computed_set_name,
    ):
        self.user_id = user_id
        self.sdf_filename = sdf_filename
        self.submitter_name = submitter_name
        self.submitter_method = submitter_method
        self.target = target
        self.version = version
        self.zfile = zfile
        self.zfile_hashvals = zfile_hashvals
        self.computed_set_name = computed_set_name

    def process_pdb(self, pdb_code, zfile, zfile_hashvals) -> str | None:
        for key in zfile_hashvals.keys():
            if key == pdb_code:
                pdb_code = f'{pdb_code}#{zfile_hashvals[pdb_code]}'

        try:
            pdb_fp = zfile[pdb_code]
        except KeyError:
            return None

        # ensure filename uniqueness
        pdb_fn = '_'.join([zfile[pdb_code].split('/')[-1], uuid.uuid4().hex])
        pdb_field = Path(settings.COMPUTED_SET_MEDIA_DIRECTORY).joinpath(pdb_fn)

        new_filename = Path(settings.MEDIA_ROOT).joinpath(pdb_field)
        old_filename = Path(settings.MEDIA_ROOT).joinpath(pdb_fp)
        old_filename.rename(new_filename)
        os.chmod(new_filename, 0o755)

        return str(pdb_field)

    # use zfile object for pdb files uploaded in zip
    def get_site_observation(
        self, property_name, mol, target, compound_set, zfile, zfile_hashvals
    ) -> SiteObservation | str | None:
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

        if zfile:
            # pdb archive uploaded. referenced pdb file may or may not be included
            pdb_code = pdb_fn.replace('.pdb', '')
            pdb_file = self.process_pdb(
                pdb_code=pdb_code,
                zfile=zfile,
                zfile_hashvals=zfile_hashvals,
            )
            if pdb_file:
                return pdb_file
            else:
                logger.info(
                    'No protein pdb (%s) found in zipfile',
                    pdb_fn,
                )

        # pdb was not included, try to find the matching site observation
        name = pdb_fn
        site_obvs = None
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

    def create_mol(self, inchi, target, name=None) -> Compound:
        # check for an existing compound, returning a Compound

        sanitized_mol = Chem.MolFromInchi(inchi, sanitize=True)
        Chem.RemoveStereochemistry(sanitized_mol)
        inchi = Chem.inchi.MolToInchi(sanitized_mol)
        inchi_key = Chem.InchiToInchiKey(inchi)

        try:
            # NB! Max said there could be thousands of compounds per
            # target so this distinct() here may become a problem

            # fmt: off
            cpd = Compound.objects.filter(
                computedmolecule__computed_set__target=target,
            ).distinct().get(
                inchi_key=inchi_key,
            )
            # fmt: on
        except Compound.DoesNotExist:
            cpd = Compound(
                smiles=Chem.MolToSmiles(sanitized_mol),
                inchi=inchi,
                inchi_key=inchi_key,
                current_identifier=name,
            )
            cpd.save()
        except MultipleObjectsReturned as exc:
            # NB! when processing new uploads, Compound is always
            # fetched by inchi_key, so this shouldn't ever create
            # duplicates. Ands LHS uploads do not create inchi_keys,
            # so under normal operations duplicates should never
            # occur. However there's nothing in the db to prevent
            # this, so adding a catch clause and writing a meaningful
            # message
            logger.error(
                'Duplicate compounds for target %s with inchi key %s.',
                target.title,
                inchi_key,
            )
            raise MultipleObjectsReturned from exc

        return cpd

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

        compound: Compound = self.create_mol(
            inchi, compound_set.target, name=molecule_name
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

                ref = qs.order_by('-cmpd_id').first()

            insp_frags.append(ref)

        ref_property = 'ref_pdb'
        # data in ref ref_pdb field may be one of 2 things:
        # - siteobservation's short code (code field)
        # - pdb file in uploaded zipfile
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

        # Need a ComputedMolecule before saving.
        # Check if anything exists already...

        # I think, realistically, I only need to check compound
        # fmt: off
        qs = ComputedMolecule.objects.filter(
            compound=compound,
        ).annotate(
            # names come in format:
            # target_name-sequential number-sequential letter,
            # e.g. A71EV2A-1-a, hence grabbing the 3rd column
            suffix=Func(
                F('name'),
                Value('-'),
                Value(3),
                function='split_part',
                output_field=TextField(),
            ),
        )

        if qs.exists():
            suffix = next(
                alphanumerator(start_from=qs.order_by('-suffix').first().suffix)
            )
        else:
            suffix = 'a'

        # distinct is ran on indexed field, so shouldn't be a problem
        number = ComputedMolecule.objects.filter(
            computed_set__target=compound_set.target,
        ).values('id').distinct().count() + 1
        # fmt: on

        name = f'v{number}{suffix}'

        existing_computed_molecules = []
        for k in qs:
            kmol = Chem.MolFromMolBlock(k.sdf_info)
            if kmol:
                # find distances between corresponding atoms of the
                # two conformers. if any one exceeds the _DIST_LIMIT,
                # consider it to be a new ComputedMolecule
                _, _, atom_map = Chem.rdMolAlign.GetBestAlignmentTransform(mol, kmol)
                molconf = mol.GetConformer()
                kmolconf = kmol.GetConformer()
                small_enough = True
                for mol_atom, kmol_atom in atom_map:
                    molpos = np.array(molconf.GetAtomPosition(mol_atom))
                    kmolpos = np.array(kmolconf.GetAtomPosition(kmol_atom))
                    distance = np.linalg.norm(molpos - kmolpos)
                    if distance >= _DIST_LIMIT:
                        small_enough = False
                        break
                if small_enough:
                    existing_computed_molecules.append(k)

        if len(existing_computed_molecules) == 1:
            logger.warning(
                'Using existing ComputedMolecule %s and overwriting its metadata',
                existing_computed_molecules[0],
            )
            computed_molecule = existing_computed_molecules[0]
        elif len(existing_computed_molecules) > 1:
            logger.warning('Deleting existing ComputedMolecules (more than 1 found')
            for exist in existing_computed_molecules:
                logger.info('Deleting ComputedMolecule %s', exist)
                exist.delete()
            computed_molecule = ComputedMolecule(name=name)
        else:
            logger.info('Creating new ComputedMolecule')
            computed_molecule = ComputedMolecule(name=name)

        if isinstance(ref_so, SiteObservation):
            code = ref_so.code
            pdb_info = ref_so.experiment.pdb_info
            lhs_so = ref_so
        else:
            code = None
            pdb_info = ref_so
            lhs_so = None

        # I don't quite understand why the overwrite of existing
        # compmol.. but this is how it was, not touching it now
        # update: I think it's about updating metadata. moving
        # name attribute out so it won't get overwritten
        computed_molecule.compound = compound
        computed_molecule.sdf_info = Chem.MolToMolBlock(mol)
        computed_molecule.site_observation_code = code
        computed_molecule.reference_code = code
        computed_molecule.molecule_name = molecule_name
        computed_molecule.smiles = smiles
        computed_molecule.pdb = lhs_so
        # TODO: this is wrong
        computed_molecule.pdb_info = pdb_info
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

        compound_set.computed_molecules.add(computed_molecule)

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
        try:
            computed_set = ComputedSet.objects.get(name=self.computed_set_name)
            # refresh some attributes
            computed_set.md_ordinal = F('md_ordinal') + 1
            computed_set.upload_date = datetime.date.today()
            computed_set.save()
        except ComputedSet.DoesNotExist:
            # no, create new

            today: datetime.date = datetime.date.today()
            new_ordinal: int = 1
            try:
                target = Target.objects.get(title=self.target)
            except Target.DoesNotExist as exc:
                # probably wrong target name supplied
                logger.error('Target %s does not exist', self.target)
                raise Target.DoesNotExist from exc

            cs_name: str = (
                f'{truncated_submitter_method}-{str(today)}-'
                + f'{get_column_letter(new_ordinal)}'
            )
            logger.info('Creating new ComputedSet "%s"', cs_name)

            computed_set = ComputedSet(
                name=cs_name,
                md_ordinal=new_ordinal,
                upload_date=today,
                method=self.submitter_method[: ComputedSet.LENGTH_METHOD],
                target=target,
                spec_version=float(self.version.strip('ver_')),
            )
            if self.user_id:
                try:
                    computed_set.owner_user = User.objects.get(id=self.user_id)
                except User.DoesNotExist as exc:
                    logger.error('User %s does not exist', self.user_id)
                    raise User.DoesNotExist from exc

            else:
                # The User ID may only be None if AUTHENTICATE_UPLOAD is False.
                # Here the ComputedSet owner will take on a default (anonymous) value.
                assert settings.AUTHENTICATE_UPLOAD is False

            computed_set.save()

        # check compound set folder exists.
        cmp_set_folder = os.path.join(
            settings.MEDIA_ROOT, settings.COMPUTED_SET_MEDIA_DIRECTORY
        )
        if not os.path.isdir(cmp_set_folder):
            logger.info('Making ComputedSet folder (%s)', cmp_set_folder)
            os.mkdir(cmp_set_folder)

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
