from dataclasses import dataclass

from typing import Any
# from typing import Optional
from typing import Iterable
from typing import Dict

import logging
from pathlib import Path
import hashlib
import yaml
import tarfile
import uuid
from tempfile import TemporaryDirectory

import functools

from celery import Task

# from hypothesis.definitions import VectTypes

# from hypothesis.models import Vector3D
# from hypothesis.models import Vector

# from frag.network.decorate import get_3d_vects_for_mol

from django.core.exceptions import MultipleObjectsReturned

from django.utils import timezone

from django.db import IntegrityError
from django.db import transaction


from django.contrib.auth import get_user_model

from scoring.models import SiteObservationGroup

from viewer.models import (
    Compound,
    Project,
    Xtalform,
    XtalformSite,
    QuatAssembly,
    CanonSite,
    CanonSiteConf,
    SiteObservation,
    Experiment,
    ExperimentUpload,
    XtalformQuatAssembly,
    Target,
    SiteObservationTag,
    TagCategory,
)

from fragalysis.settings import TARGET_LOADER_MEDIA_DIRECTORY

from django.conf import settings

logger = logging.getLogger(__name__)

# data that goes to tables are in the following files
METADATA_FILE = "meta_aligner.yaml"
ASSEMBLIES_FILE = "assemblies.yaml"
XTALFORMS_FILE = "xtalforms.yaml"
ASSIGNED_XTALFORMS_FILE = "assigned_xtalforms.yaml"
CONFIG_FILE = "config*.yaml"


# data blocks from from meta_aligner.yaml are processed into dictionaries:
# { some_id: MetadataObjects, ...}
# reason being, quite often I need to refer to these by
# some alternative id. with the dataclass, I'm able
# to create temporary dicts with key that's needed atm.
# 'new' as a flag for newly created objects (as opposed to fetched
# from the database)
@dataclass
class MetadataObjects:
    instance: Any
    data: Any = None
    new: bool = True


def create_object(func):
    """Decorator function for saving database objects"""

    @functools.wraps(func)
    def wrapper_create_object(*args, **kwargs):
        logger.debug("+ wrapper_service_query")
        logger.debug("args passed: %s", args)
        logger.debug("kwargs passed: %s", kwargs)
        logger.debug("function: %s", func.__name__)

        try:
            obj, new = func(*args, **kwargs)
        except IntegrityError as exc:
            # msg = f"Failed to save Experiment: {protein_name}"
            # update_task(self.task, "ERROR", msg)
            # logger.error("%s%s", self.task_id, msg)
            raise IntegrityError from exc
        except ValueError as exc:
            # update_task(self.task, "ERROR", exc.args[0])
            raise ValueError(exc.args[0]) from exc

        # return MetadataObjects(
        #     instance=obj, data=data.get("aligned_files", None), new=new
        # )

        return
    return wrapper_create_object


class TargetLoader:
    def __init__(
        self,
        data_bundle: str,
        proposal_ref: str,
        tempdir: str,
        user_id=None,
        task=None,
    ):
        self.data_bundle = Path(data_bundle).name
        self.bundle_name = Path(data_bundle).stem
        self.bundle_path = data_bundle
        self.proposal_ref = proposal_ref
        self.tempdir = tempdir
        self.raw_data = Path(self.tempdir).joinpath(self.bundle_name)
        self.task = task
        self.version_number = None
        self.version_dir = None
        self.previous_version_dirs = None
        self.user_id = user_id

        self.raw_data.mkdir()

        # create exp upload object
        # NB! this is not saved here in case upload fails
        self.experiment_upload = ExperimentUpload(
            commit_datetime=timezone.now(),
            file=self.data_bundle,
        )

        # work out where the data finally lands
        # path = Path(settings.MEDIA_ROOT).joinpath(TARGET_LOADER_DATA)
        path = Path(TARGET_LOADER_MEDIA_DIRECTORY)

        # give each upload a unique directory. since I already have
        # task_id, why not reuse it
        if task:
            path = path.joinpath(str(task.request.id))
            self.experiment_upload.task_id = task.request.id
        else:
            # unless of course I don't have task..
            # TODO: i suspect this will never be used.
            path_uuid = uuid.uuid4().hex
            path = path.joinpath(path_uuid)
            self.experiment_upload.task_id = path_uuid

        # figure out absolute and relative paths to final
        # location. relative path is added to db field, this will be
        # used in url requests to retrieve the file. absolute path is
        # for moving the file to the final location
        self._final_path = path.joinpath(self.bundle_name)
        self._abs_final_path = (
            Path(settings.MEDIA_ROOT).joinpath(path).joinpath(self.bundle_name)
        )
        # but don't create now, this comes later

        # to be used in logging messages, if no task, means invoked
        # directly, likely from management command
        self.task_id = f"task {task.request.id}: " if task else ""

        # these will be filled later
        self.target_name = None
        self._target_root = None
        self.target = None
        self.project = None

    @property
    def target_root(self) -> Path:
        return self._target_root

    @property
    def final_path(self) -> Path:
        return self._final_path

    @property
    def abs_final_path(self) -> Path:
        return self._abs_final_path

    def _load_yaml(self, yaml_file: Path) -> dict:
        try:
            with open(yaml_file, "r", encoding="utf-8") as file:
                contents = yaml.safe_load(file)
        except FileNotFoundError as exc:
            msg = f"{yaml_file.stem} file not found in data archive!"
            logger.error("%s%s", self.task_id, msg)
            raise FileNotFoundError(msg) from exc

        return contents

    def _get_meta_blocks(self, yaml_file, blocks) -> list:
        validation_errors = []
        error_text = "'{}' section missing in input file"
        result = []
        for block in blocks:
            try:
                result.append(yaml_file[block])
            except KeyError:
                msg = error_text.format(block)
                validation_errors.append(msg)
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)

        # wonder if it's worth throwing a custom exception here.. easier
        # to scan the logs. then again, if errors are passed to user,
        # inspecting logs isn't really necessary?
        if validation_errors:
            raise KeyError(validation_errors)

        return result

    @create_object
    def create_experiment(self, protein_name=None, data=None):
        files = self._check_file_struct(
            self.target_root, data["crystallographic_files"]
        )

        # what I would want is smth like this:
        pdb_info, mtz_info, cif_info = self.validate_files(
            file_struct=data["crystallographic_files"],
            required=("xtal_pdb", "mtz_pdb"),
            recommended=("cif_info",),
        )

        # for f in ("xtal_pdb", "mtz_pdb"):
        #     validate_file

        pdb_info = self._get_final_path(files["xtal_pdb"])
        mtz_info = self._get_final_path(files["mtz_pdb"])

        # TODO:this may be missing from the struct
        cif_info = self._get_final_path(files.get("cif_info", None))

        exp_type = (1 if data["type"] == "manual" else 0,)  # FIXME

        # status	int	new/deprecated/superseded	metadata.yaml
        # version	int	old versions are kept	target loader

        status = 1
        version = 1

        return Experiment.objects.get_or_create(
            experiment_upload=self.experiment_upload,
            code=protein_name,
            defaults={
                "status": status,
                "version": version,
                "type": exp_type,
                "pdb_info": pdb_info,
                "mtz_info": mtz_info,
                "cif_info": cif_info,
                # this doesn't seem to be present
                # pdb_sha256:
            },
        )

    def _process_experiment(self, existing_objects=None, protein_name=None, data=None):
        """Create Experiment model instance from data.

        Incoming data format:
        type: <type>
        last_updated: <datetime>
        crystallographic_files:
          xtal_pdb: {file: <file path>, sha256: <hash>}
          xtal_mtz: {file: <file path>, sha256: <hash>}
          ligand_cif: {file: <file path>, sha256: <hash>, smiles: <smiles>}
          panddas_event_files:
          - {file: <file path>, sha256: <hash>,
            model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
          - {file: <file path>, sha256: <hash>,
            model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        status: <status>
        assigned_xtalform: <str>
        aligned_files: [...]

        Manages to save all references to other tables
        """
        logger.debug("Creating experiment object: %s", protein_name)

        try:
            experiment = existing_objects.get(code=protein_name)
            new = False
        except Experiment.DoesNotExist:
            new = True
            files = self._check_file_struct(
                self.target_root, data["crystallographic_files"]
            )
            experiment = Experiment(
                experiment_upload=self.experiment_upload,
                code=protein_name,
                status=1,
                version=1,
                type=1 if data["type"] == "manual" else 0,  # FIXME
                pdb_info=str(self._get_final_path(files["xtal_pdb"])),
                mtz_info=str(self._get_final_path(files["xtal_mtz"])),
                # this may be missing from the struct
                cif_info=str(self._get_final_path(files.get("cif_info", None))),
                # this doesn't seem to be present
                # pdb_sha256=
            )
            try:
                experiment.save()
            except IntegrityError as exc:
                msg = f"Failed to save Experiment: {protein_name}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc
            except ValueError as exc:
                update_task(self.task, "ERROR", exc.args[0])
                raise ValueError(exc.args[0]) from exc

        return MetadataObjects(
            instance=experiment, data=data.get("aligned_files", None), new=new
        )

    def _process_compound(self, existing_objects=None, protein_name=None, data=None):
        """Create Compound model instance from data.

        Incoming data format:
        xtal_pdb: {file: <file path>, sha256: <hash>}
        xtal_mtz: {file: <file path>, sha256: <hash>}
        ligand_cif: {file: <file path>, sha256: <hash>, smiles: <smiles>}
        panddas_event_files:
        - {file: <file path>, sha256: <hash>,
          model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        - {file: <file path>, sha256: <hash>,
          model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        """

        # TODO: there's a method, calc_cpd, could I use that? as I
        # understand, they don't need to be unique anymore, so I can
        # just remove the uniqueness validation?

        smiles = data["ligand_cif"]["smiles"]
        try:
            compound = existing_objects.get(smiles=smiles)
            new = False
        except Compound.DoesNotExist:
            new = True

            compound = Compound(smiles=smiles)

            try:
                compound.save()
            except IntegrityError as exc:
                msg = f"Failed to save Compound: {protein_name}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc
            except KeyError as exc:
                # this means ligand info missing
                raise KeyError from exc

            compound.project_id.add(self.experiment_upload.project)

        # data basically just contains file locations, need to copy them
        # later
        return MetadataObjects(instance=compound, data=data, new=new)

    def _process_xtalform(self, existing_objects=None, idx=None, data=None):
        """Create Xtalform model instance from data.

        Incoming data format (from meta_aligner.yaml):
        <name>:
          xtalform_ref: <ref>
          xtalform_space_group: <space group>
          xtalform_cell: <cell info>

        and (from xtalforms.yaml):
        <name>:
            reference: <ref>
            assemblies:
                <idx>:
                    assembly: <assembly_id>
                    chains: <chains>

        Saves all references to other tables (QuatAssembly and Experiment).
        """
        logger.debug("Creating Xtalform object: %s", data["xtalform_ref"])
        try:
            xtalform = existing_objects.get(name=idx)
            new = False
        except Xtalform.DoesNotExist:
            new = True
            xtalform = Xtalform(
                name=idx,
                space_group=data["xtalform_space_group"],
                unit_cell_info=data["xtalform_cell"],
            )
            try:
                xtalform.save()
            except IntegrityError as exc:
                msg = f"Failed to save Xtalform: {data['xtalform_ref']}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

        return MetadataObjects(instance=xtalform, data=None, new=new)

    def _process_quat_assembly(self, existing_objects=None, idx=None, data=None):
        """Create QuatAssemblylform model instance from data.

        Incoming data format:
        <idx>:
            reference: <name>
            biomol: <biomol: str>
            chains: <chain info: str>

        No references to other models.
        """
        logger.debug("Creating QuatAssembly object: %d", idx)
        try:
            quat_assembly = existing_objects.get(chains=data["chains"])
            new = False
        except QuatAssembly.DoesNotExist:
            new = True

            quat_assembly = QuatAssembly(
                chains=data["chains"],
                name=idx,
            )
            try:
                quat_assembly.save()
            except IntegrityError as exc:
                msg = f"Failed to save QuatAssembly: {idx}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=quat_assembly,
            data=None,
            new=new,
        )

    def _process_canon_site_conf(
        self, existing_objects=None, canon_site=None, idx=None, data=None
    ):
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx: str>:
          reference_ligand_id: <lig_ref>
          residues: <array[char]>
          members: <array[char]>

        Unable to add references to:
        - SiteObservation (ref_site_observation)
        """
        logger.debug("Creating CanonSiteConf object: %s", idx)
        try:
            canon_site_conf = existing_objects.get(name=idx)
            new = False
        except CanonSiteConf.DoesNotExist:
            new = True

            canon_site_conf = CanonSiteConf(
                name=idx,
                residues=data["residues"],
                canon_site=canon_site,
            )

            try:
                canon_site_conf.save()
            except IntegrityError as exc:
                msg = f"Failed to save CanonSiteConf: {data['name']}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=canon_site_conf,
            data=data["reference_ligand_id"],
            new=new,
        )

    def _process_site_observation(
        self,
        existing_objects=None,
        experiment=None,
        compound=None,
        xtalform_site=None,
        canon_site_conf=None,
        chain=None,
        ligand=None,
        idx=None,
        data=None,
    ):
        """Create SiteObservation model instance from data.

        Incoming data format:
        <idx (apparently canon_site_conf id)>: {
          structure: <file path>,
          artefacts:  <file path>,
          event_map:  <file path>,
          x_map:  <file path>,
          pdb_apo:  <file path>,
          pdb_apo_solv:  <file path>,
          pdb_apo_desolv:  <file path>,
          ligand_mol:  <file path>,
          ligand_pdb:  <file path>,
          ligand_smiles: <smiles>,
        }
        """
        code = f"{experiment.code}_{chain}_{str(ligand)}_{str(idx)}"
        logger.debug("Creating SiteObservation object: %s", code)
        try:
            site_observation = existing_objects.get(code=code)
            new = False
        except SiteObservation.DoesNotExist:
            new = True

            files = self._check_file_struct(self.target_root, data)

            mol_data = None
            with open(
                self.target_root.joinpath(files["ligand_mol"]),
                "r",
                encoding="utf-8",
            ) as f:
                mol_data = f.read()

            try:
                site_observation = SiteObservation(
                    # Code for this protein (e.g. Mpro_Nterm-x0029_A_501_0)
                    code=code,
                    experiment=experiment,
                    cmpd=compound,
                    xtalform_site=xtalform_site,
                    canon_site_conf=canon_site_conf,
                    bound_file=str(self._get_final_path(files["structure"])),
                    apo_solv_file=str(self._get_final_path(files["pdb_apo_solv"])),
                    apo_desolv_file=str(self._get_final_path(files["pdb_apo_desolv"])),
                    apo_file=str(self._get_final_path(files["pdb_apo"])),
                    xmap_2fofc_file=str(self._get_final_path(files["2Fo-Fc_map"])),
                    xmap_fofc_file=str(self._get_final_path(files["Fo-Fc_map"])),
                    event_file=str(self._get_final_path(files["event_map"])),
                    # artefacts file currently missing, hence the get
                    artefacts_file=str(
                        self._get_final_path(files.get("artefacts", None))
                    ),
                    pdb_header_file="currently missing",
                    smiles=data["ligand_smiles"],
                    seq_id=ligand,
                    chain_id=chain,
                    ligand_mol_file=mol_data,
                )
            except KeyError as exc:
                logger.debug("exc: %s", exc)
                msg = (
                    f"Reference to {exc} file missing from aligned_files/"
                    + f"{experiment.code}/{chain}/{str(ligand)}/{str(idx)}"
                )
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise KeyError(msg) from exc

            try:
                site_observation.save()
            except IntegrityError as exc:
                msg = f"Failed to save SiteObservation: {site_observation.code}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

            # TODO: may have to implement KeyError, in case file is
            # missing. but that's only when it's guaranteed that it
            # should exist.

        return MetadataObjects(
            instance=site_observation,
            data=None,
            new=new,
        )

    def _process_canon_site(self, existing_objects=None, idx=None, data=None):
        """Create CanonSite model instance from data.

        Incoming data format:
        <id: str>:
            conformer_site_ids: <array[str]>
            global_reference_dtag: <str>
            reference_conformer_site_id: <str>
            residues: <array[str]>

        Unable to add references to:
        - CanonSiteConf (ref_conf_site)

        """
        logger.debug("Creating CanonSite object: %d", idx)
        try:
            canon_site = existing_objects.get(name=idx)
            new = False
        except CanonSite.DoesNotExist:
            new = True

            canon_site = CanonSite(
                name=idx,
                residues=data["residues"],
            )
            try:
                canon_site.save()
            except IntegrityError as exc:
                msg = f"Failed to save CanonSite: {idx}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=canon_site,
            data=data,
            new=new,
        )

    def _process_xtalform_quatassembly(
        self,
        existing_objects=None,
        xtalform=None,
        quat_assembly=None,
        idx=None,
        data=None,
    ):
        """Create XtalformQuatAssembly model instance from data.

        Incoming data format:
        <idx>:
            assembly: <assembly id: int>
            chains: <str>

        """
        logger.debug("Creating XtalformQuatAssembly object: %d", idx)
        try:
            xtal_quat = existing_objects.get(
                xtalform=xtalform,
                quat_assembly=quat_assembly,
                assembly_id=idx,
            )
            new = False
        except XtalformQuatAssembly.DoesNotExist:
            new = True

            xtal_quat = XtalformQuatAssembly(
                xtalform=xtalform,
                quat_assembly=quat_assembly,
                chains=data["chains"],
                assembly_id=idx,
            )
            try:
                xtal_quat.save()
            except IntegrityError as exc:
                msg = f"Failed to save XtalformQuatAssembly: {idx}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=xtal_quat,
            data=data,
            new=new,
        )

    def _process_xtalform_site(
        self, existing_objects=None, xtalform=None, canon_site=None, idx=None, data=None
    ):
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx>:
          xtalform_id: <str>
          canonical_site_id: <str>
          crystallographic_chain: A
          members: <array[str]>

        Saves references to all other tables (Xtalform and CanonSite).
        """
        logger.debug("Creating XtalformSite object: %d", idx)
        try:
            xtalform_site = existing_objects.get(xtalform_site_id=idx)
            new = False
        except XtalformSite.DoesNotExist:
            new = True

            xtalform_site = XtalformSite(
                xtalform=xtalform,
                canon_site=canon_site,
                lig_chain=data["crystallographic_chain"],
                residues=data["members"],
                xtalform_site_id=idx,
            )

            try:
                xtalform_site.save()
            except IntegrityError as exc:
                msg = f"Failed to save Xtalform: {idx}"
                update_task(self.task, "ERROR", msg)
                logger.error("%s%s", self.task_id, msg)
                raise IntegrityError(msg) from exc

        return MetadataObjects(
            instance=xtalform_site,
            data=data["members"],
            new=new,
        )

    def process_metadata(
        self,
        upload_root: Path = None,
        task: Task = None,
    ):
        """Extract model instances from metadata file and save them to db."""
        # TODO: this method is quite long and should perhaps be broken
        # apart. then again, it just keeps calling other methods to
        # create model instances, so logically it's just doing the
        # same thing.

        update_task(task, "PROCESSING", "Processing metadata")
        logger.info("%sProcessing %s", self.task_id, upload_root)

        # moved this bit from init
        self.target, target_created = Target.objects.get_or_create(
            title=self.target_name
        )

        # TODO: original target loader's function get_create_projects
        # seems to handle more cases. adopt or copy
        visit = self.proposal_ref.split()[0]
        self.project, project_created = Project.objects.get_or_create(title=visit)

        try:
            committer = get_user_model().objects.get(pk=self.user_id)
        except get_user_model().DoesNotExist:
            # add upload as anonymous user
            committer = get_user_model().objects.get(pk=settings.ANONYMOUS_USER)

        # TODO: is it here where I can figure out if this has already been uploaded?
        if self._is_already_uploaded(target_created, project_created):
            # remove uploaded file
            Path(self.bundle_path).unlink()
            raise FileExistsError(f"{self.bundle_name} already uploaded, skipping.")

        if project_created and committer.pk == settings.ANONYMOUS_USER:
            self.project.open_to_public = True
            self.project.save()

        # populate m2m field
        self.target.project_id.add(self.project)

        self.experiment_upload.project = self.project
        self.experiment_upload.target = self.target
        self.experiment_upload.committer = committer
        self.experiment_upload.save()

        assemblies = self._load_yaml(Path(upload_root).joinpath(ASSEMBLIES_FILE))
        xtalform_assemblies = self._load_yaml(
            Path(upload_root).joinpath(XTALFORMS_FILE)
        )
        meta = self._load_yaml(Path(upload_root).joinpath(METADATA_FILE))
        assigned_xtalforms = self._load_yaml(
            Path(upload_root).joinpath(ASSIGNED_XTALFORMS_FILE)
        )

        # collect top level info
        self.version_number = meta["version_number"]
        self.version_dir = meta["version_dir"]
        self.previous_version_dirs = meta["previous_version_dirs"]

        blocks = [
            "crystals",
            "xtalforms",
            "canon_sites",
            "conformer_sites",
            "xtalform_sites",
        ]
        try:
            (  # pylint: disable=unbalanced-tuple-unpacking
                crystals,
                xtalforms,
                canon_sites,
                conformer_sites,
                xtalform_sites,
            ) = self._get_meta_blocks(meta, blocks)
        except FileNotFoundError as exc:
            raise FileNotFoundError(exc.args[0]) from exc
        except KeyError as exc:
            raise KeyError(exc.args[0]) from exc

        result = []

        # memo to self: the order of saving objects is dictated by db
        # relations - when handling version 2 upload, I need to
        # occasionally check if the object already exists and if yes,
        # fetch it from the db. cannot do this without quering target
        # (or visit?) so need to save objects along the relationships

        # fetch existing objects
        old_experiments = Experiment.filter_manager.by_target(self.target)
        logger.debug("%s existing Experiment objects found", old_experiments.count())

        old_compounds = Compound.filter_manager.by_target(self.target)
        logger.debug("%s existing Compound objects found", old_compounds.count())

        old_xtalforms = Xtalform.filter_manager.by_target(self.target)
        logger.debug("%s existing Xtalform objects found", old_xtalforms.count())

        old_xtalquatasm = XtalformQuatAssembly.filter_manager.by_target(self.target)
        logger.debug(
            "%s existing XtalformQuatAssembly objects found", old_xtalquatasm.count()
        )

        old_quatassemblies = QuatAssembly.filter_manager.by_target(self.target)
        logger.debug(
            "%s existing QuatAssembly objects found", old_quatassemblies.count()
        )

        old_xtalformsites = XtalformSite.filter_manager.by_target(self.target)
        logger.debug(
            "%s existing XtalformSite objects found", old_xtalformsites.count()
        )

        old_canonsites = CanonSite.filter_manager.by_target(self.target)
        logger.debug("%s existing CanonSite objects found", old_canonsites.count())

        old_canonsiteconfs = CanonSiteConf.filter_manager.by_target(self.target)
        logger.debug(
            "%s existing CanonSiteConf objects found", old_canonsiteconfs.count()
        )

        old_siteobservations = SiteObservation.filter_manager.by_target(self.target)
        logger.debug(
            "%s existing SiteObservaton objects found", old_siteobservations.count()
        )

        experiment_objects = {}
        compound_objects = {}
        for prot_name, prot_data in crystals.items():
            # TODO: unclear if I need to save experiment that doesn't have
            # aligned crystal files section or should I move this beyond
            # continue as well.
            try:
                experiment_objects[prot_name] = self._process_experiment(
                    existing_objects=old_experiments,
                    protein_name=prot_name,
                    data=prot_data,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

            cmpd_data = prot_data["crystallographic_files"]
            try:
                compound_objects[prot_name] = self._process_compound(
                    existing_objects=old_compounds,
                    protein_name=prot_name,
                    data=cmpd_data,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc
            except KeyError:
                # this particular block doesn't have compound info
                # continue with the loop, nothing to do
                continue

        result.append(self._log_msg(experiment_objects))
        result.append(self._log_msg(compound_objects))

        # save components manytomany to experiment
        # TODO: is it 1:1 relationship? looking at the meta_align it seems to be,
        # but why the m2m then?
        for comp_code, comp_meta in compound_objects.items():
            experiment = experiment_objects[comp_code].instance
            experiment.compounds.add(comp_meta.instance)

        xtalform_objects = {}
        for idx, obj in xtalforms.items():
            try:
                xtalform_objects[idx] = self._process_xtalform(
                    existing_objects=old_xtalforms,
                    idx=idx,
                    data=obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

        result.append(self._log_msg(xtalform_objects))

        # add xtalform fk to experiment
        for _, obj in experiment_objects.items():
            try:
                obj.instance.xtalform = xtalform_objects[
                    assigned_xtalforms[obj.instance.code]
                ].instance
                obj.instance.save()
            except KeyError:
                # TODO: don't know at this point whether ignoring this
                # error is a right thing to do
                msg = f"xtalform {obj.instance.code} undefined for {obj}"
                result.append(msg)
                logger.warning(msg)

        quat_assembly_objects = {}
        for idx, obj in assemblies.items():
            try:
                quat_assembly_objects[idx] = self._process_quat_assembly(
                    existing_objects=old_quatassemblies,
                    idx=idx,
                    data=obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

        result.append(self._log_msg(quat_assembly_objects))

        # this is used just for logging, no other function
        xtalform_quat_assembly_objects = {}
        for xtalform_id, data in xtalform_assemblies.items():
            xtalform = xtalform_objects[xtalform_id].instance
            for idx, obj in data["assemblies"].items():
                quat_assembly = quat_assembly_objects[obj["assembly"]].instance
                key = f"{xtalform_id} {quat_assembly.id} {idx}"
                xtalform_quat_assembly_objects[
                    key
                ] = self._process_xtalform_quatassembly(
                    existing_objects=old_xtalquatasm,
                    xtalform=xtalform,
                    quat_assembly=quat_assembly,
                    idx=idx,
                    data=obj,
                )

        result.append(self._log_msg(xtalform_quat_assembly_objects))

        canon_site_objects = {}
        for idx, obj in canon_sites.items():
            try:
                canon_site_objects[idx] = self._process_canon_site(
                    existing_objects=old_canonsites,
                    idx=idx,
                    data=obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc
        # NB! missing fk's:
        # - ref_conf_site
        # - quat_assembly

        result.append(self._log_msg(canon_site_objects))

        # reindex canon sites by canon_sites_conf_sites
        # NB! this is also used below for ref_conf_site in canon_site
        canon_sites_by_conf_sites = {
            conf: obj.instance
            for obj in canon_site_objects.values()
            for conf in obj.data["conformer_site_ids"]
        }

        canon_site_conf_objects = {}
        for idx, obj in conformer_sites.items():
            try:
                canon_site_conf_objects[idx] = self._process_canon_site_conf(
                    existing_objects=old_canonsiteconfs,
                    canon_site=canon_sites_by_conf_sites[idx],
                    idx=idx,
                    data=obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc
        # NB! missing fk's:
        # - site_observation

        result.append(self._log_msg(canon_site_conf_objects))

        xtalform_sites_objects = {}
        for idx, obj in xtalform_sites.items():
            try:
                xtalform_sites_objects[idx] = self._process_xtalform_site(
                    existing_objects=old_xtalformsites,
                    xtalform=xtalform_objects[obj["xtalform_id"]].instance,
                    canon_site=canon_site_objects[obj["canonical_site_id"]].instance,
                    idx=idx,
                    data=obj,
                )
            except IntegrityError as exc:
                raise IntegrityError(exc.args[0]) from exc

        result.append(self._log_msg(xtalform_sites_objects))

        # now can update CanonSite with ref_conf_site
        for val in canon_site_objects.values():
            val.instance.ref_conf_site = canon_site_conf_objects[
                val.data["reference_conformer_site_id"]
            ].instance
            val.instance.save()

        # canon site instances are now complete
        # still missing fk to site_observation in canon_site_conf

        # reindex xtalform site to grab for site observation
        xtalform_site_by_tag = {}
        for val in xtalform_sites_objects.values():
            for k in val.data:
                xtalform_site_by_tag[k] = val.instance

        site_observation_objects = {}
        # TODO: would be nice to get rid of quadruple for
        for experiment_meta in experiment_objects.values():
            if experiment_meta.data is None:
                continue
            for chain, ligand in experiment_meta.data.items():
                for ligand, ligand_data in ligand.items():
                    for idx, obj in ligand_data.items():
                        key = f"{experiment_meta.instance.code}/{chain}/{ligand}"
                        try:
                            site_observation_objects[
                                key
                            ] = self._process_site_observation(
                                existing_objects=old_siteobservations,
                                experiment=experiment_meta.instance,
                                compound=compound_objects[
                                    experiment_meta.instance.code
                                ].instance,
                                xtalform_site=xtalform_site_by_tag[key],
                                canon_site_conf=canon_site_conf_objects[idx].instance,
                                chain=chain,
                                ligand=ligand,
                                idx=idx,
                                data=obj,
                            )
                        except IntegrityError as exc:
                            raise IntegrityError(exc.args[0]) from exc

        result.append(self._log_msg(site_observation_objects))

        tag_categories = (
            "ConformerSites",
            "CanonSites",
            "XtalformSites",
            "Quatassemblies",
            "Xtalforms",
        )

        for cat in tag_categories:
            self._tag_site_observations(site_observation_objects, cat)

        # final remaining fk, attach reference site observation to canon_site_conf
        for val in canon_site_conf_objects.values():
            val.instance.ref_site_observation = site_observation_objects[
                val.data
            ].instance
            val.instance.save()

        result.append(
            f"{self.bundle_name} {upload_root.name}: "
            + f"User {self.experiment_upload.committer} "
            + f"uploaded target {self.target}"
        )

        return result

    def _tag_site_observations(self, site_observation_objects, category):
        # this is an attempt to replicate tag creation from previous
        # loader. as there are plenty of differences, I cannot just
        # use the same functions..

        logger.debug("Getting category %s", category)
        groups = {}
        for _, obj in site_observation_objects.items():
            if category == "ConformerSites":
                tags = [
                    "conf_site: " + ",".join(obj.instance.canon_site_conf.residues[:3]),
                ]
            elif category == "CanonSites":
                tags = [
                    "canon_site: "
                    + ",".join(obj.instance.xtalform_site.canon_site.residues[:3]),
                ]
            elif category == "XtalformSites":
                tags = [
                    "xtalform_site: "
                    + ",".join(obj.instance.xtalform_site.residues[:3]),
                ]
            # tricky one. connected via m2m
            elif category == "Quatassemblies":
                tags = [
                    "quatassembly: " + qa.name
                    for qa in obj.instance.xtalform_site.xtalform.quat_assembly.all()
                ]

            elif category == "Xtalforms":
                tags = [
                    "xtalform: " + obj.instance.xtalform_site.xtalform.name,
                ]
            else:
                tags = [
                    "Unspecified",
                ]

            for tag in tags:
                if tag not in groups.keys():
                    groups[tag] = [obj.instance]
                else:
                    groups[tag].append(obj.instance)

        # I suspect I need to group them by site..
        for tag, so_list in groups.items():
            try:
                # memo to self: description is set to tag, but there's
                # no fk to tag, instead, tag has a fk to
                # group. There's no uniqueness requirement on
                # description so there's no certainty that this will
                # be unique (or remain searchable at all because user
                # is allowed to change the tag name). this feels like
                # poor design but I don't understand the principles of
                # this system to know if that's indeed the case or if
                # it is in fact a truly elegant solution
                so_group = SiteObservationGroup.objects.get(
                    target=self.target, description=tag
                )
            except SiteObservationGroup.DoesNotExist:
                so_group = SiteObservationGroup(target_id=self.target.pk)
                so_group.save()
            except MultipleObjectsReturned:
                SiteObservationGroup.objects.filter(
                    target=self.target, description=tag
                ).delete()
                so_group = SiteObservationGroup(target_id=self.target.pk)
                so_group.save()

            try:
                so_tag = SiteObservationTag.objects.get(tag=tag, target=self.target)
                # Tag already exists
                # Apart from the new mol_group and molecules, we shouldn't be
                # changing anything.
                so_tag.mol_group = so_group
            except SiteObservationTag.DoesNotExist:
                so_tag = SiteObservationTag()
                so_tag.tag = tag
                so_tag.category = TagCategory.objects.get(category=category)
                so_tag.target = self.target
                so_tag.mol_group = so_group

            so_tag.save()

            for site_obvs in so_list:
                logger.debug("site_obvs_id=%s", site_obvs.id)
                so_group.site_observation.add(site_obvs)
                so_tag.site_observations.add(site_obvs)

    def _is_already_uploaded(self, target_created, project_created):
        if target_created or project_created:
            return False
        else:
            uploaded_files = ExperimentUpload.objects.filter(
                target=self.target,
                project=self.project,
            ).values_list("file", flat=True)

            # TODO: this just tests the target-project-filename combo,
            # which may not be enough

            return self.data_bundle in uploaded_files

    def _get_final_path(self, path: Path):
        """Update relative path to final storage path

        NB! this returns a relative path that can be used in queries
        not absoulte one. This is used to populate location fields in
        database tables.
        """
        try:
            return self.final_path.joinpath(self.target_name).joinpath(path)
        except TypeError:
            # received invalid path
            return None

    def process_bundle(self, task=None):
        """Resolves subdirs in uploaded data bundle.

        If called from task, takes task as a parameter for status updates.
        """
        # result is final report, list of messages passed back to user
        result = []

        # by now I should have archive unpacked, get target name from
        # config.yaml

        # this a bit of an chicken and an egg problem. I need to get
        # to the config file in upload_N directory, but the path to
        # that goes through another, which is usually the same as
        # target. But.. I don't know if I can trust that, which is why
        # I grab it from the config file where says this is the target
        it = self.raw_data.iterdir()
        try:
            target_dir = next(it)
        except StopIteration as exc:
            raise StopIteration("Target directory missing from uploaded file!") from exc

        logger.debug("target_dir: %s", target_dir)

        # a quick sanity check, assuming only one target per upload
        if sum(1 for _ in it) != 0:
            tgd = "; ".join([str(p) for p in it])
            raise AssertionError(
                f"More than one target directory in uploaded file: {tgd}"
            )

        target_path = self.raw_data.joinpath(target_dir)
        it = target_path.iterdir()
        try:
            upload_dir = next(it)
        except StopIteration as exc:
            raise StopIteration("Upload directory missing from uploaded file!") from exc

        # technically, what I could do here is the same validation as
        # with targed dir above, there should be only on upload_N
        # directory under target. but it doesn't matter if there's
        # more, I'll just use the first one

        config_it = upload_dir.glob(CONFIG_FILE)
        try:
            config_file = next(config_it)
        except StopIteration as exc:
            raise StopIteration(f"config file missing from {upload_dir}") from exc

        config = self._load_yaml(config_file)
        logger.debug("config: %s", config)

        try:
            self.target_name = config["target_name"]
        except KeyError as exc:
            raise KeyError("target_name missing in config file!") from exc

        self._target_root = self.raw_data.joinpath(self.target_name)

        # as mentioned above, there should be only one upload_N
        # directory, so I suppose I could get rid of the loop
        # here. doesn't hurt though, so I'll leave it for now
        for path in Path(self.target_root).iterdir():
            if path.is_dir():
                logger.info("Found upload dir: %s", str(path))
                try:
                    upload_report = self.process_metadata(
                        upload_root=path,
                        task=task,
                    )
                except FileNotFoundError as exc:
                    result.append(exc.args[0])
                    raise FileNotFoundError(exc.args[0]) from exc
                except IntegrityError as exc:
                    result.append(exc.args[0])
                    raise IntegrityError(exc.args[0]) from exc
                except ValueError as exc:
                    result.append(exc.args[0])
                    raise ValueError(exc.args[0]) from exc
                except FileExistsError as exc:
                    # this target has been uploaded at- least once and
                    # this upload already exists. skip
                    result.append(exc.args[0])
                    raise FileExistsError(exc.args[0]) from exc

                result.extend(upload_report)

        return result

    # standardized logging when processing metadata file
    def _log_msg(self, obj_dict):
        new_obj = sum([1 for k in obj_dict.values() if k.new])
        item = next(iter(obj_dict.values()))
        msg = f"{len(obj_dict.keys())} {item.instance._meta.model} objects processed, {new_obj} created"  # pylint: disable=protected-access
        update_task(self.task, "PROCESSING", msg)
        logger.info("%s%s", self.task_id, msg)
        return f"{self.bundle_name} {self.version_dir}: {msg}"

    def validate_files(
        self,
        file_struct: Dict,
        required: Iterable[str] = (),
        recommended: Iterable[str] = (),
    ) -> list[str]:
        """Check if file exists and if sha256 hash matches (if given).

        file struct can come in 2 configurations:
        {file_key: {file: <file_path>, sha265: <hash> [smiles: <smiles>]}, ...}
        or simply
        {file_key: <file path>}
        Detect which one and take appropriate action.

        F
        """
        result = {}
        for key, value in file_struct.items():
            filename, file_hash = None, None

            # sort out the filename
            if isinstance(value, dict):
                file_hash = value.get("sha256", None)
                try:
                    filename = value["file"]
                except KeyError as exc:
                    # this is rather unexpected, haven't seen it yet

                    # NB! I only have the key at this point so I
                    # cannot reference it properly in here. I think I
                    # need to raise a custom exception here
                    msg = "Malformed dict, key 'file' missing!"
                    if key in required:
                        logger.error(msg)
                        raise KeyError(msg) from exc
                    elif key in recommended:
                        logger.warning(msg)
                        continue
            elif isinstance(value, str) and not key == "ligand_smiles":
                # workaround for ligand_smiles, it happens to be
                # specified on the same level as files
                filename = value
            else:
                # this is probably the list of panddas event files, don't
                # need them here
                # although.. should i validate them here nevertheless?
                # i'd have to do this on copy otherwise..
                continue

            # now I have file path and possibly hash
            file_path = self.target_root.joinpath(filename)
            if file_path.is_file():
                if file_hash:
                    if file_hash == self._calculate_sha256(file_path):
                        result[key] = filename
                    else:
                        msg =f"Invalid hash for file {filename}!"
                        if key in required:
                            logger.error(msg)
                            raise ValueError(msg)
                        elif key in recommended:
                            logger.warning(msg)
            else:
                msg = f"{key} referenced in {METADATA_FILE} but not found in archive!"
                if key in required:
                    logger.error(msg)
                    raise FileNotFoundError(msg)
                elif key in recommended:
                    logger.warning(msg)

        files = []
        # for f in required:
        #     try:
        #         files.append(result[f])
        #     except KeyError as exc
        #         raise KeyError from exc

        # for f in recommended:
        #     try:
        #         files.append(result[f])
        #     except KeyError:
        #         raise KeyError

        return [ result[k] for k in required ] + [ result[k] for k in recommended ]

    # trying to get rid of that
    @staticmethod
    def _check_file_struct(upload_root, file_struct):
        """Check if file exists and if sha256 hash matches (if given).

        file struct can come in 2 configurations:
        {file_key: {file: <file_path>, sha265: <hash> [smiles: <smiles>]}, ...}
        or simply
        {file_key: <file path>}
        Detect which one and take appropriate action.
        """
        result = {}
        # for key, value in file_struct.items():
        #     if isinstance(value, dict):
        #         try:
        #             filename = value["file"]

        #             # I guess I don't care if it's not given?
        #             file_hash = value.get("sha256", None)

        #             # TODO: error handling. don't know how to
        #             # implement this because data structure/metadata
        #             # not at final structure yet
        #             if TargetLoader._check_file(
        #                 upload_root.joinpath(filename), file_hash
        #             ):
        #                 result[key] = str(filename)

        #         except KeyError:
        #             logger.warning("%s file info missing in %s!", key, METADATA_FILE)
        #         except ValueError as exc:
        #             logger.error("Invalid hash for file %s!", filename)
        #             raise ValueError(f"Invalid hash for file {filename}!") from exc

        #     elif isinstance(value, str) and not key == "ligand_smiles":
        #         # this is a bit of a workaround but ligand_smiles
        #         # happens to be specified on the same level as files
        #         if TargetLoader._check_file(upload_root.joinpath(value)):
        #             result[key] = str(value)
        #         else:
        #             # wait, I think I should actually raise exp here..
        #             logger.warning(
        #                 "%s referenced in %s but not found in archive",
        #                 value,
        #                 METADATA_FILE,
        #             )
        #     else:
        #         # this is probably the list of panddas event files, don't
        #         # need them here
        #         # although.. should i validate them here nevertheless?
        #         # i'd have to do this on copy otherwise..
        #         pass

        return result

    # @staticmethod
    # def _check_file(file_path: Path, file_hash=None):
    #     """Check if file exist and compare with hash."""
    #     if file_path.is_file():
    #         if file_hash:
    #             if file_hash == TargetLoader._calculate_sha256(file_path):
    #                 return True
    #             else:
    #                 # not logging error here because don't have
    #                 # correct file path
    #                 raise ValueError
    #         else:
    #             return True
    #     return False


    # borrowed from SO
    @staticmethod
    def _calculate_sha256(filepath):
        sha256_hash = hashlib.sha256()
        with open(filepath, "rb") as f:
            # Read the file in chunks of 4096 bytes
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()


def load_target(
    data_bundle,
    proposal_ref=None,
    contact_email=None,
    user_id=None,
    task=None,
):
    # TODO: do I need to sniff out correct archive format?
    with TemporaryDirectory(dir=settings.MEDIA_ROOT) as tempdir:
        target_loader = TargetLoader(
            data_bundle, proposal_ref, tempdir, user_id=user_id, task=task
        )
        try:
            # archive is first extracted to temporary dir and moved later
            with tarfile.open(target_loader.bundle_path, "r") as archive:
                msg = f"Extracting bundle: {data_bundle}"
                logger.info("%s%s", target_loader.task_id, msg)
                update_task(task, "PROCESSING", msg)
                archive.extractall(target_loader.raw_data)
                msg = f"Data extraction complete: {data_bundle}"
                logger.info("%s%s", target_loader.task_id, msg)
        except FileNotFoundError as exc:
            msg = f"{data_bundle} file does not exist!"
            logger.exception("%s%s", target_loader.task_id, msg)
            target_loader.experiment_upload.message = exc.args[0]
            raise FileNotFoundError(msg) from exc

        try:
            with transaction.atomic():
                upload_report = target_loader.process_bundle(task=task)
        except FileNotFoundError as exc:
            logger.error(exc.args[0])
            target_loader.experiment_upload.message = exc.args[0]
            raise FileNotFoundError(exc.args[0]) from exc
        except IntegrityError as exc:
            logger.error(exc, exc_info=True)
            target_loader.experiment_upload.message = exc.args[0]
            raise IntegrityError(exc.args[0]) from exc
        except ValueError as exc:
            logger.error(exc, exc_info=True)
            raise IntegrityError(exc.args[0]) from exc
        except FileExistsError as exc:
            logger.error(exc.args[0])
            target_loader.experiment_upload.message = exc.args[0]
            raise FileExistsError(exc.args[0]) from exc
        except AssertionError as exc:
            logger.error(exc.args[0])
            target_loader.experiment_upload.message = exc.args[0]
            raise AssertionError(exc.args[0]) from exc

        # move to final location
        target_loader.abs_final_path.mkdir(parents=True)
        target_loader.raw_data.rename(target_loader.abs_final_path)
        Path(target_loader.bundle_path).rename(
            target_loader.abs_final_path.parent.joinpath(target_loader.data_bundle)
        )

        update_task(task, "SUCCESS", upload_report)
        target_loader.experiment_upload.message = upload_report

        logger.debug("%s", upload_report)

        target_loader.experiment_upload.save()


def update_task(task, state, message):
    try:
        task.update_state(
            state=state,
            meta={
                "description": message,
            },
        )
    except AttributeError:
        # no task passed to method, nothing to do
        pass
