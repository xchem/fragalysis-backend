import ast
import copy
import datetime
import logging
import os
import re
import uuid
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from dateutil.parser import parse
from openpyxl.utils import get_column_letter

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
import django
from django.db import IntegrityError, transaction

django.setup()

from django.conf import settings
from django.core.exceptions import MultipleObjectsReturned, ValidationError
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.validators import validate_email
from django.db.models import F
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

from .sdf_check import add_warning

logger = logging.getLogger(__name__)


# maximum distance between corresponding atoms in poses
_DIST_LIMIT = 0.5

EMPTY_VALUES = (
    "nan",
    "",
    None,
    np.nan,
)


HEADER_MOL_FIELDS = (
    "ref_url",
    "method",
    "submitter_name",
    "submitter_institution",
    "submitter_email",
    "generation_date",
)


# How do we get the 'prefix' and 'version' from the MOL Name.
# (used at the moment to handle mol_refs in Squonk-generated SD files).
# They look like this: -
#   A71EV2A-x0379_A_147_1_A71EV2A-x0379+A+147+1_LIG
#   ------------------- -
#      "Prefix"      "Version"
# And we want a 'long code' from this, e.g.: -
#   A71EV2A-x0379_A_147_v1
_re_ref_mol_long_code = re.compile(
    r"(?P<prefix>([^_]*)_(\S+)_(\d+))_(?P<version>\d+)_(.*)"
)


def dataType(a_str: str) -> str:
    lean_str = a_str.strip()
    if not lean_str:
        return "BLANK"

    try:
        t = ast.literal_eval(lean_str)
    except (ValueError, SyntaxError):
        return "TEXT"
    else:
        if type(t) in [int, int, float, bool]:
            if t in [
                True,
                False,
                "TRUE",
                "FALSE",
                "true",
                "false",
                "yes",
                "no",
                "YES",
                "NO",
                "Yes",
                "No",
                "Y",
                "N",
                "y",
                "n",
            ]:
                return "BIT"
            if type(t) is int or type(t) is int:
                return "INT"
            if type(t) is float:
                return "FLOAT"

            # Can't get here?
            assert False
        else:
            return "TEXT"


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
                if filename.split(".")[-1] == "pdb":
                    # Test if Protein object already exists
                    code = filename.split("/")[-1].replace(".pdb", "")
                    test_pdb_code = filename.split("/")[-1].replace(".pdb", "")
                    test_site_obvs_objs = SiteObservation.objects.filter(
                        code=test_pdb_code
                    )

                    if len(test_site_obvs_objs) != 0:
                        # make a unique pdb code as not to overwrite existing object
                        rand_str = uuid.uuid4().hex
                        test_pdb_code = f"{code}#{rand_str}"
                        zfile_hashvals[code] = rand_str

                    fn = f"{test_pdb_code}.pdb"
                    pdb_path = default_storage.save(
                        f"tmp/{fn}", ContentFile(zf.read(filename))
                    )
                    zfile[test_pdb_code] = pdb_path

            # Close the zip file
            if zf:
                zf.close()

        return zfile, zfile_hashvals

    def run(self, params) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, str]]]:
        return self.save_pdb_zip(params["pdb_zip"])


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
        computed_set_id,
    ):
        self.user_id = user_id
        self.sdf_filename = sdf_filename
        self.submitter_name = submitter_name
        self.submitter_method = submitter_method
        self.target_id = target
        self.version = version
        self.zfile = zfile
        self.zfile_hashvals = zfile_hashvals
        self.computed_set_id = computed_set_id

        # using the same mechanism to pass messages as validation
        self.messages: dict[str, Any] = {
            "molecule_name": [],
            "field": [],
            "warning_string": [],
        }

    def process_pdb(self, pdb_code, zfile, zfile_hashvals) -> str | None:
        for key in zfile_hashvals.keys():
            if key == pdb_code:
                pdb_code = f"{pdb_code}#{zfile_hashvals[pdb_code]}"

        try:
            pdb_fp = zfile[pdb_code]
        except KeyError:
            return None

        # ensure filename uniqueness
        pdb_fn = "_".join([zfile[pdb_code].split("/")[-1], uuid.uuid4().hex])
        pdb_field = Path(settings.COMPUTED_SET_MEDIA_DIRECTORY).joinpath(pdb_fn)

        new_filename = Path(settings.MEDIA_ROOT).joinpath(pdb_field)
        old_filename = Path(settings.MEDIA_ROOT).joinpath(pdb_fp)

        # there may be a case where 2 or more molfiles reference the
        # same pdb. in this case, the old pdb is already renamed to
        # new.
        if old_filename.exists() and not new_filename.exists():
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

        pdb_fn = mol.GetProp(property_name).split("/")[-1]

        if zfile:
            # pdb archive uploaded. referenced pdb file may or may not be included
            pdb_code = pdb_fn.replace(".pdb", "")
            pdb_file = self.process_pdb(
                pdb_code=pdb_code,
                zfile=zfile,
                zfile_hashvals=zfile_hashvals,
            )
            if pdb_file:
                return pdb_file
            else:
                logger.info(
                    "No protein pdb (%s) found in zipfile",
                    pdb_fn,
                )

        # pdb was not included, try to find the matching site observation
        name = pdb_fn
        site_obvs = None
        try:
            site_obvs = SiteObservation.objects.get(
                code__contains=name,
                experiment__experiment_upload__target__pk=target,
                superseded=False,
            )
        except SiteObservation.DoesNotExist:
            # Initial SiteObservation lookup failed.
            logger.warning(
                "Failed to get SiteObservation object (target=%s name=%s)",
                compound_set.target.title,
                name,
            )
            # Try alternatives.
            # If all else fails then the site_obvs will be 'None'
            qs = SiteObservation.objects.filter(
                code__contains=name,
                experiment__experiment_upload__target__pk=target,
            )
            if qs.exists():
                logger.info(
                    "Found SiteObservation containing name=%s qs=%s",
                    name,
                    qs,
                )
            else:
                alt_name = name.split(":")[0].split("_")[0]
                qs = SiteObservation.objects.filter(
                    code__contains=alt_name,
                    experiment__experiment_upload__target__pk=target,
                )
                if qs.exists():
                    logger.info(
                        "Found SiteObservation containing alternative name=%s qs=%s",
                        alt_name,
                        qs,
                    )
            if qs.count() > 0:
                logger.debug(
                    "Found alternative (target=%s name=%s)",
                    compound_set.target.title,
                    name,
                )
                site_obvs = qs[0]

        if not site_obvs:
            logger.warning(
                "No SiteObservation found (target=%s pdb_fn=%s)",
                compound_set.target.title,
                pdb_fn,
            )

        return site_obvs

    def create_mol(self, inchi, target, name=None) -> tuple[Compound, str]:
        # check for an existing compound, returning a Compound

        sanitized_mol = Chem.MolFromInchi(inchi, sanitize=True)
        Chem.RemoveStereochemistry(sanitized_mol)
        inchi = Chem.inchi.MolToInchi(sanitized_mol)
        inchi_key = Chem.InchiToInchiKey(inchi)

        # look for *all* compounds under this target, both LHS and RHS
        # uploads
        lhs_qs = Compound.filter_manager.by_target(target)
        rhs_qs = Compound.objects.filter(
            pk__in=ComputedMolecule.objects.filter(
                computed_set__target=target,
            ).values('compound')
        )

        cpd_number = "1"
        try:
            cpd = rhs_qs.get(inchi_key=inchi_key)
            # memo to self: I'm not setting cpd_number here, because
            # it's read from computedmol name
        except Compound.DoesNotExist:
            # RHS didn't work out, try LHS. here, duplicates are possible
            cpd = lhs_qs.filter(inchi_key=inchi_key).first()

            if not cpd:
                # still no compound, create new
                cpd = Compound(
                    smiles=Chem.MolToSmiles(sanitized_mol),
                    inchi=inchi,
                    inchi_key=inchi_key,
                    description=name,
                )
                # This is a new compound.
                cpd.save()
                # This is a new compound.
                # We must now set relationships to the Proposal that it applies to.
                cpd.project_id.add(target.project)
                qs = Compound.objects.filter(
                    computedmolecule__computed_set__target=target,
                )
                cpd_number = str(qs.count())
        except MultipleObjectsReturned as exc:
            # NB! when processing new uploads, Compound is always
            # fetched by inchi_key, so this shouldn't ever create
            # duplicates. Ands LHS uploads do not create inchi_keys,
            # so under normal operations duplicates should never
            # occur. However there's nothing in the db to prevent
            # this, so adding a catch clause and writing a meaningful
            # message
            msg = f"Duplicate compounds for target {target.title} with inchi key {inchi_key}."
            logger.error(msg)
            raise IntegrityError(msg) from exc

        return cpd, cpd_number

    def set_props(self, cpd, props, score_descriptions):
        for sd, val in score_descriptions.items():
            logger.debug("sd: %s", sd)
            logger.debug("sd.name, val: %s: %s", sd.name, val)
            if dataType(str(props[sd.name])) == "TEXT":
                score_value = TextScoreValues()
            else:
                score_value = NumericalScoreValues()

            try:
                float(val)
            except ValueError:
                return None

            if sd.name in HEADER_MOL_FIELDS:
                score_value.value = val
            else:
                score_value.value = props[sd.name]

            score_value.compound = cpd
            score_value.score = sd
            score_value.save()

        return None

    def set_mol(
        self, mol, target, compound_set, filename, zfile=None, zfile_hashvals=None
    ) -> ComputedMolecule:
        # Don't need...
        assert target
        assert compound_set

        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.inchi.MolToInchi(mol)
        molecule_name = mol.GetProp("_Name")

        flattened_copy = copy.deepcopy(mol)
        Chem.RemoveStereochemistry(mol)
        flat_inchi = Chem.inchi.MolToInchiKey(flattened_copy)
        logger.debug('flattened inchi key: %s', flat_inchi)

        compound, number = self.create_mol(
            inchi, compound_set.target, name=molecule_name
        )

        insp = mol.GetProp("ref_mols")
        insp = insp.split(",")
        insp = [i.strip() for i in insp]
        insp_frags = []
        for i in insp:
            # try exact match first
            logger.debug('looking for so code %s', str(i))
            try:
                site_obvs = SiteObservation.objects.get(
                    code=str(i),
                    experiment__experiment_upload__target=compound_set.target,
                    superseded=False,
                )
                ref = site_obvs
            except SiteObservation.DoesNotExist:
                # A hack - for Squonk Job execution tests.
                # The ref_mols field doesn't contain a long or short code to simplify lookups.
                # To get to a long code we can use the defined reg-ex pattern.
                long_code = ""
                re_match = _re_ref_mol_long_code.match(i)
                if re_match:
                    prefix = re_match.group('prefix')
                    version_number = re_match.group('version')
                    # Long code is the 'prefix' and the 'version' (with a 'v')
                    long_code = f"{prefix}_v{version_number}"
                if not long_code:
                    raise IntegrityError(  # pylint: disable=raise-missing-from
                        f"Could not find long-code pattern in {i}"
                    )
                logger.warning(
                    "Search for '%s' failed - now looking for '%s' (target=%s)...",
                    i,
                    long_code,
                    compound_set.target,
                )
                try:
                    qs = SiteObservation.objects.filter(
                        longcode=long_code,
                        experiment__experiment_upload__target=compound_set.target,
                    )
                except SiteObservation.DoesNotExist:
                    raise IntegrityError(  # pylint: disable=raise-missing-from
                        f"No matching molecules found for inspiration frag {i}"
                    )

                # Why order-by here and not above?
                # And - is the response ever not 0 or 1 records?
                ref = qs.order_by("-cmpd_id").first()

            insp_frags.append(ref)

        ref_property = "ref_pdb"
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
                "Failed to get a Reference SiteObservation (%s) for %s, %s, %s",
                ref_property,
                mol,
                target,
                compound_set,
            )

        # Need a ComputedMolecule before saving.
        # Check if anything exists already...

        # I think, realistically, I only need to check compound
        # update: I used to annotate name components, with the new
        # format, this is not necessary. or possible
        qs = ComputedMolecule.objects.filter(
            compound=compound,
        ).order_by("name")

        if qs.exists():
            # not actually latest, just last according to sorting above
            latest = qs.last()
            # regex pattern - split name like 'v1a'
            # ('(letters)(digits)(letters)' to components
            groups = re.search(r"()(\d+)(\D+)", qs.last().name)
            if groups is None or len(groups.groups()) != 3:
                # just a quick sanity check
                raise IntegrityError(
                    f"Non-standard ComputedMolecule.name: {latest.name}"
                )
            number = groups.groups()[1]  # type: ignore [index]
            suffix = next(alphanumerator(start_from=groups.groups()[2]))  # type: ignore [index]
        else:
            suffix = "a"
            # number = 1

        name = f"v{number}{suffix}"

        existing_computed_molecules = []
        for k in qs:
            if k.compound.inchi_key == flat_inchi:
                # existing compound is a flattened copy of the new one, match found
                existing_computed_molecules.append(k)
                continue

            kmol = Chem.MolFromMolBlock(k.sdf_info)
            if kmol:
                # find distances between corresponding atoms of the
                # two conformers. if any one exceeds the _DIST_LIMIT,
                # consider it to be a new ComputedMolecule
                try:
                    _, _, atom_map = Chem.rdMolAlign.GetBestAlignmentTransform(
                        mol, kmol
                    )
                except RuntimeError as exc:
                    msg = (
                        f"Failed to find alignment between {k.molecule_name} "
                        + f'and {mol.GetProp("original ID")}'
                    )
                    logger.error(msg)
                    raise IntegrityError(msg) from exc

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
                "Using existing ComputedMolecule %s and overwriting its metadata",
                existing_computed_molecules[0],
            )
            computed_molecule = existing_computed_molecules[0]
        elif len(existing_computed_molecules) > 1:
            logger.warning("Deleting existing ComputedMolecules (more than 1 found")
            for exist in existing_computed_molecules:
                logger.info("Deleting ComputedMolecule %s", exist)
                exist.delete()
            computed_molecule = ComputedMolecule(name=name)
        else:
            logger.info("Creating new ComputedMolecule (name=%s)", name)
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
        # compmol ... but this is how it was, not touching it now
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
            mol.GetProp("ref_url") if mol.HasProp("ref_url") else None
        )
        computed_molecule.ref_url = ref_url if is_url(ref_url) else None
        rationale: Optional[str] = (
            mol.GetProp("rationale") if mol.HasProp("rationale") else None
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
        datestring = description_mol.GetProp("generation_date")
        try:
            date = parse(datestring, dayfirst=True)
        except ValueError as exc:
            msg = f'"{datestring}" is not a valid date'
            logger.error(msg)
            raise IntegrityError(msg) from exc

        submitter, created = ComputedSetSubmitter.objects.get_or_create(
            method=description_mol.GetProp("method"),
            email=description_mol.GetProp("submitter_email"),
        )
        if created:
            submitter.name = description_mol.GetProp("submitter_name")
            submitter.institution = description_mol.GetProp("submitter_institution")
            submitter.generation_date = date
            submitter.save()

        return submitter

    def process_mol(
        self,
        mol,
        target,
        compound_set,
        filename,
        score_descriptions,
        zfile=None,
        zfile_hashvals=None,
    ) -> None:
        molecule_name = mol.GetProp("_Name")
        logger.debug("+ process_mol %s", molecule_name)

        other_props = mol.GetPropsAsDict()
        skip_mol = False

        # if ref_mols or ref_pdb is missing skip the molecule
        for prop in ["ref_mols", "ref_pdb"]:
            if prop not in other_props.keys():
                self.messages = add_warning(
                    molecule_name=molecule_name,
                    field=prop,
                    warning_string=f"Property {prop} missing. Skipping molecule!",
                    validate_dict=self.messages,
                )
                skip_mol = True
            elif other_props[prop] in EMPTY_VALUES:
                self.messages = add_warning(
                    molecule_name=molecule_name,
                    field=prop,
                    warning_string=f"Property {prop} undefined. Skipping molecule!",
                    validate_dict=self.messages,
                )
                skip_mol = True

        # if any header mol fields are defined on non-header molecules those values are ignored and a warning shown
        for prop in HEADER_MOL_FIELDS:
            if prop not in other_props.keys():
                # non-header molecules don't need header fields
                continue

            if other_props[prop] not in EMPTY_VALUES:
                # header fields in non-header molecules have values ignored
                self.messages = add_warning(
                    molecule_name=molecule_name,
                    field=prop,
                    warning_string=f"Property {prop} value {other_props[prop]} ignored.",
                    validate_dict=self.messages,
                )

            # get rid of the header field property on the non-header molecule
            del other_props[prop]

        if skip_mol:
            logger.warning("Skipping molecule '%s'", molecule_name)
        else:
            cpd = self.set_mol(
                mol, target, compound_set, filename, zfile, zfile_hashvals
            )
            self.set_props(cpd, other_props, score_descriptions)

    def set_descriptions(
        self, filename, computed_set: ComputedSet
    ) -> tuple[List[Chem.rdchem.Mol], dict[str, ScoreDescription]]:
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
        if description_mol.HasProp("ref_url"):
            computed_set.method_url = description_mol.GetProp("ref_url")
        computed_set.save()

        description_dict = description_mol.GetPropsAsDict()
        logger.debug("index mol original values: %s", description_dict)
        # score descriptions for this upload, doesn't matter if
        # created or existing
        score_descriptions = {}
        errors = []
        for key in description_dict.keys():
            if key in descriptions_needed and key not in [
                "ref_mols",
                "ref_pdb",
                "index",
                "Name",
            ]:
                description, _ = ScoreDescription.objects.get_or_create(
                    computed_set=computed_set,
                    name=key,
                    description=description_dict[key],
                )

                value = description_dict[key]

                if key in HEADER_MOL_FIELDS:
                    if value in EMPTY_VALUES:
                        msg = f"Empty value for {key} in header molecule"
                        errors.append(msg)
                        logger.error(msg)
                    if key == "submitter_email":
                        try:
                            validate_email(value)
                        except ValidationError:
                            msg = f'"{value}" is not a valid email'
                            logger.error(msg)
                            errors.append(msg)

                score_descriptions[description] = value

        logger.debug("index mol values: %s", score_descriptions.values())
        if errors:
            raise IntegrityError(errors)

        return mols, score_descriptions

    def task(self) -> tuple[ComputedSet, dict]:
        # Truncate submitted method (lower-case)?
        truncated_submitter_method: str = "unspecified"
        try:
            with transaction.atomic():
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
                    computed_set = ComputedSet.objects.get(pk=self.computed_set_id)
                    # refresh some attributes
                    computed_set.md_ordinal = F("md_ordinal") + 1
                    computed_set.upload_date = datetime.date.today()
                    computed_set.save()
                except (ValueError, ComputedSet.DoesNotExist):
                    # ValueError when pk is None
                    # no, create new

                    today: datetime.date = datetime.date.today()
                    new_ordinal: int = 1

                    try:
                        target = Target.objects.get(pk=self.target_id)
                    except Target.DoesNotExist as exc:
                        # target's existance should be validated in the view,
                        # this could hardly happen
                        msg = f"Target {self.target_id} does not exist"
                        logger.error(msg)
                        raise IntegrityError(msg) from exc

                    cs_name: str = (
                        f"{truncated_submitter_method}-{str(today)}-"
                        + f"{get_column_letter(new_ordinal)}"
                    )

                    # now that I have a name, I can check whether this
                    # target already has this set
                    try:
                        # this feels wrong, I think it's better if the
                        # object is resolved in the view.. or maybe in
                        # validate task..
                        computed_set = ComputedSet.objects.get(
                            name=cs_name,
                            target=target,
                        )
                    except ComputedSet.DoesNotExist:
                        # and only now create new set
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
                            msg = f"User {self.user_id} does not exist"
                            logger.error(msg)
                            raise IntegrityError(msg) from exc

                    else:
                        # The User ID may only be None if AUTHENTICATE_UPLOAD is False.
                        # Here the ComputedSet owner will take on a default (anonymous) value.
                        assert settings.AUTHENTICATE_UPLOAD is False

                    computed_set.save()

                # Set descriptions in return for the Molecules.
                # This also sets the submitter and method URL properties of the computed set
                # while also saving it.
                sdf_filename = str(self.sdf_filename)
                mols_to_process, score_descriptions = self.set_descriptions(
                    filename=sdf_filename, computed_set=computed_set
                )

                # Process the molecules
                logger.info("%s mols_to_process=%s", computed_set, len(mols_to_process))
                for i in range(len(mols_to_process)):
                    logger.debug(
                        "processing mol %s: %s", i, mols_to_process[i].GetProp("_Name")
                    )
                    self.process_mol(
                        mols_to_process[i],
                        self.target_id,
                        computed_set,
                        sdf_filename,
                        score_descriptions,
                        self.zfile,
                        self.zfile_hashvals,
                    )
        except IntegrityError as exc:
            # clean up previously written files. this is not ideal,
            # they should be written to a tempdir or something, like
            # in target loader. TODO for later
            try:
                for p in self.zfile.values():
                    Path(p).unlink()
            except AttributeError:
                # zfile is None, nothing to do
                pass

            raise ValueError(exc.args[0]) from exc

        # assuming no errors, write the files

        # check compound set folder exists.
        cmp_set_folder = os.path.join(
            settings.MEDIA_ROOT, settings.COMPUTED_SET_MEDIA_DIRECTORY
        )
        if not os.path.isdir(cmp_set_folder):
            logger.info("Making ComputedSet folder (%s)", cmp_set_folder)
            os.mkdir(cmp_set_folder)

        # move and save the compound set
        new_filename = (
            Path(settings.MEDIA_ROOT)
            .joinpath(settings.COMPUTED_SET_MEDIA_DIRECTORY)
            .joinpath(
                f"{computed_set.name}_upload_{computed_set.md_ordinal}_{Path(sdf_filename).name}"
            )
        )
        os.rename(sdf_filename, new_filename)
        computed_set.submitted_sdf = Path(sdf_filename).name
        computed_set.written_sdf_filename = new_filename
        computed_set.save()

        logger.info("Created %s", computed_set)

        return computed_set, self.messages


def blank_mol_vals(sdf_file) -> Tuple[str, str, str]:
    """Returns the submitter name, method and version (_Name) if present.
    If not present the corresponding values are empty strings.
    """
    suppl = Chem.SDMolSupplier(sdf_file)
    if not suppl:
        return "", "", ""
    # print('%d mols detected (including blank mol)' % (len(suppl),))
    blank_mol = suppl[0]
    if not blank_mol:
        return "", "", ""

    # Get submitter name/info for passing into upload to get unique name
    submitter_name = ""
    if blank_mol.HasProp("submitter_name"):
        submitter_name = blank_mol.GetProp("submitter_name")

    submitter_method = ""
    if blank_mol.HasProp("method"):
        submitter_method = blank_mol.GetProp("method")

    version = ""
    if blank_mol.HasProp("_Name"):
        version = blank_mol.GetProp("_Name")

    return submitter_name, submitter_method, version
