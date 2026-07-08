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
from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.validators import validate_email
from django.db import IntegrityError, transaction
from django.db.models import F
from django.utils import timezone
from openpyxl.utils import get_column_letter
from rdkit import Chem

from viewer.models import (
    Compound,
    ComputedInspiration,
    ComputedSet,
    ComputedSetSubmitter,
    Result,
    ResultProperty,
    ResultValueDataType,
    SiteObservation,
    Target,
    User,
)
from viewer.utils import (
    add_props_to_sdf_molecule,
    alphanumerator,
    is_url,
    set_directory_permissions,
    word_count,
)

from .sdf_check import add_warning
from .tags import TagManager
from .target_loader import assign_observation_quality_status

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
    """Analyse string and return postgres-compatible data type"""
    lean_str = a_str.strip()
    if not lean_str:
        return "BLANK"

    try:
        t = ast.literal_eval(lean_str)
    except (ValueError, SyntaxError):
        return "text"
    else:
        if type(t) in [int, int, float, bool]:
            # this doesn't work, does it? or am I missing something? seems ureachable
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
                return "boolean"
            if type(t) is int or type(t) is int:
                return "INT"
            if type(t) is float:
                return "float"

            # Can't get here?
            assert False
        else:
            return "text"


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

        self.data_type_dict = {
            k.data_type: k for k in ResultValueDataType.objects.all()
        }

        self.result_properties = {}

        # This is now the third time the target is resolved. a bit
        # stupid. fix
        try:
            target = Target.objects.get(pk=self.target_id)
        except Target.DoesNotExist as exc:
            # target's existance should be validated in the view,
            # this could hardly happen
            msg = f"Target {self.target_id} does not exist"
            logger.error(msg)
            raise IntegrityError(msg) from exc

        last_upload = target.experimentupload_set.order_by('commit_datetime').last()

        # create directory for virtual files
        self.virtual_root = (
            Path(settings.TARGET_LOADER_MEDIA_DIRECTORY)
            .joinpath(
                str(target.zip_archive),
            )
            .joinpath(
                last_upload.upload_data_dir,
            )
            .joinpath(
                'virtual_files',
            )
        )

        Path(settings.MEDIA_ROOT).joinpath(self.virtual_root).mkdir(
            parents=True, exist_ok=True
        )

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

    def create_mol(self, inchi, target, name=None) -> Compound:
        # check for an existing compound, returning a Compound

        sanitized_mol = Chem.MolFromInchi(inchi, sanitize=True)
        Chem.RemoveStereochemistry(sanitized_mol)
        inchi = Chem.inchi.MolToInchi(sanitized_mol)
        inchi_key = Chem.InchiToInchiKey(inchi)

        # cpd_number = "1"
        qs = Compound.filter_manager.by_target(target)
        logger.debug('compounds by target: %s', qs.count())

        # duplicates possible
        cpd = qs.filter(inchi_key=inchi_key).first()
        logger.debug('cpd found: %s', cpd.pk if cpd else None)
        if not cpd:
            # no compound, create new
            cpd = Compound(
                smiles=Chem.MolToSmiles(sanitized_mol),
                inchi=inchi,
                inchi_key=inchi_key,
                description=name,
            )
            cpd.save()
            logger.debug('cpd not found, created new: %s', cpd.pk)
            # This is a new compound.
            # We must now set relationships to the Proposal that it applies to.
            cpd.project_id.add(target.project)

        return cpd

    def set_props(self, cpd, props, score_descriptions, computed_set) -> None:
        for property_name, val in score_descriptions.items():
            # logger.debug("score_descriptions: %s", score_descriptions)
            # logger.debug("property_name: %s", property_name)
            # logger.debug("props: %s", props)
            # logger.debug("sd.name, val: %s: %s", sd.name, val)

            data_type = self.data_type_dict.get(
                dataType(str(props[property_name])),
                self.data_type_dict['text'],
            )

            try:
                result_property = self.result_properties[property_name]
            except KeyError:
                result_property, _ = ResultProperty.objects.get_or_create(
                    result_property=property_name,
                    unit=None,
                    target=computed_set.target,
                    data_type=data_type,
                )

            if property_name in HEADER_MOL_FIELDS:
                value = val
            else:
                value = props[property_name]

            result = Result(
                raw_value=value,
                site_observation=cpd,
                result_property=result_property,
                computed_set=computed_set,
            )

            if data_type.data_type == 'float':
                result.float_value = value
            elif data_type.data_type == 'integer':
                result.int_value = value

            else:
                result.text_value = value

            result.save()

        return None

    def set_mol(
        self, mol, target, compound_set, filename, zfile=None, zfile_hashvals=None
    ) -> SiteObservation:
        # ) -> computedmolecule:
        # Don't need...
        assert target
        assert compound_set

        # the flattening part seems duplicated between create_mol
        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.inchi.MolToInchi(mol)
        molecule_name = mol.GetProp("_Name")

        flattened_copy = copy.deepcopy(mol)
        Chem.RemoveStereochemistry(mol)
        flat_inchi = Chem.inchi.MolToInchiKey(flattened_copy)
        logger.debug('flattened inchi key: %s', flat_inchi)

        # compound, number = self.create_mol(
        #     inchi, compound_set.target, name=molecule_name
        # )
        compound = self.create_mol(inchi, compound_set.target, name=molecule_name)

        insp = mol.GetProp("ref_mols")
        insp = insp.split(",")
        insp = [i.strip() for i in insp]
        logger.debug('got inspirations: %s', insp)
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

        qs = SiteObservation.objects.filter(
            cmpd=compound,
            experiment__isnull=True,  # get only virtual observations
        ).order_by("virtual_name")
        # memo to self: ordering is fine because they all should have
        # the same number component

        if qs.exists():
            # logger.debug('found existing connected compound: %s', qs)
            # not actually latest, just last according to sorting above
            latest = qs.last()
            # regex pattern - split name like 'v1a'
            # ('(letters)(digits)(letters)' to components
            groups = re.search(r"()(\d+)(\D+)", latest.virtual_name)
            if groups is None or len(groups.groups()) != 3:
                # just a quick sanity check
                raise IntegrityError(f"Non-standard virtual_name: {latest.name}")
            number = groups.groups()[1]  # type: ignore [index]
            suffix = next(alphanumerator(start_from=groups.groups()[2]))  # type: ignore [index]
        else:
            # this is getting a wrong count somehow
            # count_qs = SiteObservation.objects.filter(
            #     computed_set__isnull=False,
            # ).values(
            #     'cmpd',
            # )
            # logger.debug('did not find existing connected compound: %s', count_qs.count())
            suffix = "a"
            number = (
                SiteObservation.objects.filter(
                    computed_set__isnull=False,
                )
                .values(
                    'cmpd',
                )
                .distinct()
                .count()
                + 1
            )
            # logger.debug('but number is: %s', number)

        name = f"v{str(number)}{suffix}"

        if isinstance(ref_so, SiteObservation):
            # code = ref_so.code
            pdb_info = ref_so.experiment.pdb_info
            # lhs_so = ref_so
            xtalform_site = ref_so.xtalform_site
            canon_site_conf = ref_so.canon_site_conf
        else:
            # code = None
            pdb_info = ref_so
            # lhs_so = None
            xtalform_site = None
            canon_site_conf = None

        # computed_molecule.site_observation_code = code
        # computed_molecule.reference_code = code
        # Extract possible reference URL and Rationale
        # URLs have to be valid URLs and rationals must contain more than one word
        ref_url: Optional[str] = (
            mol.GetProp("ref_url") if mol.HasProp("ref_url") else None
        )
        ref_url = ref_url if is_url(ref_url) else None

        rationale: Optional[str] = (
            mol.GetProp("rationale") if mol.HasProp("rationale") else None
        )
        rationale = rationale if word_count(rationale) > 1 else None

        new_so = SiteObservation(
            cmpd=compound,
            xtalform_site=xtalform_site,
            canon_site_conf=canon_site_conf,
            smiles=smiles,
            virtual_name=name,
            virtual_molecule_name=molecule_name,
            virtual_ref_url=ref_url,
            virtual_rationale=rationale,
            virtual_pdb_info=pdb_info,
        )
        new_so.save()

        # identifier is auto-generated, doesn't exist until saved
        filename = self.virtual_root.joinpath(
            f"{compound_set.name}_upload_{compound_set.md_ordinal}_"
            + f"{name}_{molecule_name}_{new_so.virtual_identifier}.mol"
        )

        new_so_mol = Chem.MolToMolBlock(mol)
        sdf_filename = Path(settings.MEDIA_ROOT).joinpath(filename)
        with open(sdf_filename, "w", encoding='utf-8') as f:
            f.write(new_so_mol)

        new_so.virtual_ligand_mol = str(filename)
        new_so.save()
        assign_observation_quality_status(new_so)
        # computed_molecule.sdf_info = Chem.MolToMolBlock(mol)

        # find similar observations (former computedmolecules) and
        # add new so to the (if similar enough)

        # NB! this bit was added in 1394. the code looks like it
        # doesn't do what it was intednded to do. It did pass
        # validation though, so I'm not sure

        # existing_computed_molecules = []
        for so in qs:
            filepath = Path(settings.MEDIA_ROOT).joinpath(str(so.virtual_ligand_mol))
            so_mol = Chem.MolFromMolFile(str(filepath))
            if so_mol:
                # find distances between corresponding atoms of the
                # two conformers. if any one exceeds the _DIST_LIMIT,
                # consider it to be a new SiteObservation
                try:
                    _, _, atom_map = Chem.rdMolAlign.GetBestAlignmentTransform(
                        mol, so_mol
                    )
                except RuntimeError as exc:
                    msg = (
                        f"Failed to find alignment between {so.virtual_molecule_name} "
                        + f'and {mol.GetProp("original ID")}'
                    )
                    logger.error(msg)
                    raise IntegrityError(msg) from exc

                molconf = mol.GetConformer()
                kmolconf = so_mol.GetConformer()
                small_enough = True
                for mol_atom, so_mol_atom in atom_map:
                    molpos = np.array(molconf.GetAtomPosition(mol_atom))
                    so_molpos = np.array(kmolconf.GetAtomPosition(so_mol_atom))
                    distance = np.linalg.norm(molpos - so_molpos)
                    if distance >= _DIST_LIMIT:
                        small_enough = False
                        break
                if small_enough:
                    new_so.pose = so.pose
                    new_so.save()

        compound_set.site_observations.add(new_so)

        logger.debug('got insp_frags: %s', insp_frags)
        ComputedInspiration.objects.bulk_create(
            [
                ComputedInspiration(
                    site_observation=new_so,
                    computed_inspiration=inspiration,
                    computed_set=compound_set,
                )
                for inspiration in insp_frags
            ],
            ignore_conflicts=True,
        )

        # No update the molecule in the original file...
        add_props_to_sdf_molecule(
            sdf_file=str(sdf_filename),
            molecule=molecule_name,
            properties={"target_identifier": new_so.virtual_name},
        )

        return new_so

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
    ) -> int | None:
        molecule_name = mol.GetProp("_Name")
        logger.debug("+ process_mol %s", molecule_name)

        other_props = mol.GetPropsAsDict()
        skip_mol = False

        # logger.debug('other_props: %s', other_props)

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

        # if any header mol fields are defined on non-header molecules
        # those values are ignored and a warning shown
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
            return None
        else:
            cpd = self.set_mol(
                mol, target, compound_set, filename, zfile, zfile_hashvals
            )
            self.set_props(cpd, other_props, score_descriptions, compound_set)

            return cpd.pk

    def set_descriptions(
        self, filename, computed_set: ComputedSet
    ) -> tuple[List[Chem.rdchem.Mol], dict[str, str]]:
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
        result_descriptions = {}
        errors = []
        for key in description_dict.keys():
            if key in descriptions_needed and key not in [
                "ref_mols",
                "ref_pdb",
                "index",
                "Name",
            ]:
                # NB! change made during LHS RHS unification - don't
                # create description (result property) objects here,
                # do that in set_props when saving properties. Not the
                # best solution

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

                result_descriptions[key] = value

        logger.debug("index mol values: %s", result_descriptions.values())
        if errors:
            raise IntegrityError(errors)

        return mols, result_descriptions

    def task(self) -> tuple[ComputedSet, dict]:
        # Truncate submitted method (lower-case)?
        # truncated_submitter_method: str = "unspecified"
        try:
            with transaction.atomic():
                submitter_method: str = self.submitter_method
                if self.submitter_method:
                    submitter_method = self.submitter_method
                else:
                    submitter_method = "unspecified"
                    logger.warning(
                        'ComputedSet submitter method is not set. Using "%s"',
                        submitter_method,
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
                        f"{submitter_method}-{str(today)}-"
                        + f"{get_column_letter(new_ordinal)}"
                    )

                    # cache properties for later
                    self.result_properties = {
                        k.result_property: k
                        for k in ResultProperty.objects.filter(target=target)
                    }

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
                            # method=self.submitter_method[: ComputedSet.LENGTH_METHOD],
                            method=self.submitter_method,
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
                so_ids = []
                for i in range(len(mols_to_process)):
                    logger.debug(
                        "processing mol %s: %s", i, mols_to_process[i].GetProp("_Name")
                    )
                    so_pk = self.process_mol(
                        mols_to_process[i],
                        self.target_id,
                        computed_set,
                        sdf_filename,
                        score_descriptions,
                        self.zfile,
                        self.zfile_hashvals,
                    )
                    so_ids.append(so_pk)

                tagger = TagManager(computed_set.target, meta_category='rhs')
                so_qs = SiteObservation.objects.filter(pk__in=so_ids)
                datestr = timezone.now().date().strftime('%Y-%m-%d')
                tagger.tag_new_site_observations(
                    site_observations=so_qs,
                    new_observation_tag=f"{computed_set.name} {datestr}",
                )

                # adjust permissions for any files created
                set_directory_permissions(
                    Path(settings.MEDIA_ROOT).joinpath(self.virtual_root),
                    0o755,
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
