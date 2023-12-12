import os
import math
from dataclasses import dataclass
from dataclasses import field


from enum import Enum

from typing import Any

from typing import Optional
from typing import Iterable
from typing import TypeVar

import logging
from pathlib import Path
import hashlib
import yaml
import tarfile
import uuid
from tempfile import TemporaryDirectory

import functools

from celery import Task


from django.core.exceptions import MultipleObjectsReturned

from django.utils import timezone

from django.db import IntegrityError
from django.db import transaction
from django.db.models.base import ModelBase
from django.db.models import Model


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
# assemblies and xtalforms
XTALFORMS_FILE = "crystalforms.yaml"

# target name, nothing else
CONFIG_FILE = "config*.yaml"

# everything else
METADATA_FILE = "meta_aligner.yaml"


# type hint for Djanog model instance
ModelInstance = TypeVar("ModelInstance", bound=Model)


class UploadState(str, Enum):
    """Target loader progress state.

    PROCESSING - all good, upload in progress
    REPORTING  - upload failed, loader in reporting mode for diagnostics
    SUCCESS    - processing complete, all good
    FAILED     - processing complete, failed
    """

    PROCESSING = "PROCESSING"
    REPORTING = "REPORTING"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"


class Level(str, Enum):
    INFO = "INFO"
    WARNING = "WARNING"
    FATAL = "FATAL"


@dataclass
class MetadataObject:
    """Data structure to store freshly created model instances.

    data blocks from from meta_aligner.yaml are processed into
    dictionaries: { some_id: MetadataObjects, ...}

    Reason being, quite often I need to refer to these by some
    alternative ID. With the dataclass, I'm able to create temporary
    dicts with key that are needed.
    """

    instance: ModelInstance
    index_data: dict = field(default_factory=dict)


# type hint for wrapped yaml block processors
MetDict = TypeVar("MetDict", bound=dict[int | str, MetadataObject])


@dataclass
class ProcessedObject:
    """Data structure for creating model instances.

    Returned from methods that process yaml blocks to dictionaries
    that can be sent to django's model's get_or_create methods.
    """

    model_class: ModelBase
    fields: dict
    key: str
    defaults: dict = field(default_factory=dict)
    index_data: dict = field(default_factory=dict)
    identifier: Optional[str] = ""


@dataclass
class UploadReportEntry:
    level: Level
    message: str

    def __str__(self):
        return ": ".join([self.level, self.message])


@dataclass
class UploadReport:
    task: Task | None
    stack: list[UploadReportEntry] = field(default_factory=list)
    upload_state: UploadState = UploadState.PROCESSING
    failed: bool = False

    def __post_init__(self) -> None:
        self.task_id = f"task {self.task.request.id}: " if self.task else ""

    def log(self, level: Level, message: str) -> None:
        msg = f"{self.task_id}{message}"
        if level == Level.FATAL:
            self.failed = True
            self.upload_state = UploadState.REPORTING
            logger.error(msg)
        elif level == Level.WARNING:
            logger.warning(msg)
        else:
            # must be info
            logger.info(msg)

        self.stack.append(UploadReportEntry(level=level, message=message))
        self._update_task(message)

    def final(self, level: Level = Level.INFO, message: str = ""):
        if self.upload_state == UploadState.PROCESSING:
            self.upload_state = UploadState.SUCCESS
        else:
            self.upload_state = UploadState.FAILED

        self.stack.append(UploadReportEntry(level=level, message=message))
        self._update_task(self.json())

    def json(self):
        return [str(k) for k in self.stack]

    def _update_task(self, message: str | list):
        try:
            self.task.update_state(
                state=self.upload_state,
                meta={
                    "description": message,
                },
            )
        except AttributeError:
            # no task passed to method, nothing to do
            pass


def _flatten_dict_gen(d: dict, parent_key: tuple | str | int, depth: int):
    for k, v in d.items():
        if parent_key:
            if isinstance(parent_key, tuple):
                new_key = (*parent_key, k)
            else:
                new_key = (parent_key, k)
        else:
            new_key = k

        try:
            deep_enough = any([isinstance(x, dict) for x in v.values()])
        except AttributeError:
            continue

        if deep_enough and depth > 1:
            yield from flatten_dict(v, new_key, depth - 1)
        else:
            if isinstance(new_key, str):
                yield new_key, v
            else:
                yield *new_key, v


def flatten_dict(d: dict, parent_key: tuple | int | str = "", depth: int = 1):
    """Flatten nested dict to specified depth."""
    return _flatten_dict_gen(d, parent_key, depth)


def set_directory_permissions(path, permissions):
    for root, dirs, files in os.walk(path):
        # Set permissions for directories
        for directory in dirs:
            dir_path = os.path.join(root, directory)
            os.chmod(dir_path, permissions)

        # Set permissions for files
        for file in files:
            file_path = os.path.join(root, file)
            os.chmod(file_path, permissions)


# borrowed from SO
def calculate_sha256(filepath):
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read the file in chunks of 4096 bytes
        for chunk in iter(lambda: f.read(4096), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def create_objects(func=None, *, depth=math.inf):
    """Wrapper function for saving database objects.

    Handles common part of saving model instances, actual saving,
    logging, reporting and error handling.

    Inner functions are yaml data processing functions that return
    the model class and the data to pass to model's get_or_create
    function.

    """
    if func is None:
        return functools.partial(create_objects, depth=depth)

    @functools.wraps(func)
    def wrapper_create_objects(
        self, *args, yaml_data: dict, **kwargs
    ) -> dict[int | str, MetadataObject]:
        logger.debug("+ wrapper_service_query")
        # logger.debug("args passed: %s", args)
        logger.debug("kwargs passed: %s", kwargs)

        # flattened_data = flatten_dict(kwargs["yaml_data"], depth=depth)
        flattened_data = flatten_dict(yaml_data, depth=depth)
        result = {}
        created, existing, failed = 0, 0, 0
        for item in flattened_data:
            logger.debug("flattened data item: %s", item)
            instance_data = func(self, *args, item_data=item, **kwargs)
            logger.debug("Instance data returned: %s", instance_data)
            if not instance_data:
                continue

            obj = None
            try:
                obj, new = instance_data.model_class.filter_manager.by_target(
                    self.target
                ).get_or_create(
                    **instance_data.fields,
                    defaults=instance_data.defaults,
                )
                logger.debug(
                    "%s object %s created",
                    instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                    obj,
                )
                if new:
                    created = created + 1
                else:
                    existing = existing + 1
            except MultipleObjectsReturned:
                msg = "{}.get_or_create returned multiple objects for {}".format(
                    instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                    instance_data.fields,
                )
                self.report.log(Level.FATAL, msg)
                failed = failed + 1
            except IntegrityError:
                msg = "{} object {} failed to save".format(
                    instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                    instance_data.key,
                )
                self.report.log(Level.FATAL, msg)
                failed = failed + 1

            if not obj:
                # create fake object so I can just push the upload
                # through and compile report for user feedback
                obj = instance_data.model_class(
                    **instance_data.fields | instance_data.defaults
                )
                logger.warning(
                    "Fake %s object created: %s",
                    instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                    obj,
                )

            m = MetadataObject(instance=obj, index_data=instance_data.index_data)
            # index data here probs
            result[instance_data.key] = m

        msg = "{} {} objects processed, {} created, {} fetched from database".format(
            created + existing + failed,
            next(  # pylint: disable=protected-access
                iter(result.values())
            ).instance._meta.model,
            created,
            existing,
        )  # pylint: disable=protected-access
        self.report.log(Level.INFO, msg)

        return result

    return wrapper_create_objects


class TargetLoader:
    def __init__(
        self,
        data_bundle: str,
        proposal_ref: str,
        tempdir: str,
        user_id=None,
        task: Task | None = None,
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

        self.report = UploadReport(task=task)

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
        # self.task_id = f"task {task.request.id}: " if task else ""

        # these will be filled later
        self.target_name = None
        self._target_root = None
        self.target = None
        self.project = None

    @property
    def final_path(self) -> Path:
        return self._final_path

    @property
    def abs_final_path(self) -> Path:
        return self._abs_final_path

    def validate_files(
        self,
        obj_identifier: str,
        file_struct: dict,
        required: Iterable[str] = (),
        recommended: Iterable[str] = (),
    ) -> list[str | None]:
        """Check if file exists and if sha256 hash matches (if given).

        file struct can come in 2 configurations:
        {file_key: {file: <file_path>, sha265: <hash> [smiles: <smiles>]}, ...}
        or simply
        {file_key: <file path>}
        Detect which one and take appropriate action.

        Once the filename is extracted, check if it exists and if
        sha256 hash is given, calculate the hash and compare to the
        one in file.

        params:
        - file_struct: dictionary read from yaml file
        - required: mandatory filename keys
        - recommended: optional filename keys
        - protein_name: experiment_identifier (used for logging)

        return:
        - list of all file paths required

        Checks for 4 possible errors:
        - file is expected by the db schema but not referenced in METADATA_FILE
        - file is referenced METADATA_FILE but not present in uploaded archive
        - calculated hash doesn't match with the one in METADATA_FILE
        - dictionary in unexpected format, unable to extract filename

        """
        # TODO: don't like that I have to provide protein_name solely
        # for error logging purposes

        def logfunc(key, logstr, *args):
            if key in required:
                logger.error(logstr, *args)
                # add msg to record
                # don't i want something on class level? and maybe just override here?

            else:
                logger.warning(logstr, *args)

        result = {}
        for key, value in file_struct.items():
            if key not in required and key not in recommended:
                # schema isn't looking for this file, ignore
                continue

            filename, file_hash = None, None

            # sort out the filename
            if isinstance(value, dict):
                file_hash = value.get("sha256", None)
                try:
                    filename = value["file"]
                except KeyError:
                    # this is rather unexpected, haven't seen it yet
                    logfunc(
                        key, "%s: malformed dict, key 'file' missing!", obj_identifier
                    )

                    # unable to extract file from dict, no point to
                    # continue with hash checking
                    continue

            elif isinstance(value, str):
                filename = value
            else:
                # this is probably the list of panddas event files, don't
                # need them here
                # although.. should i validate them here nevertheless?
                # i'd have to do this on copy otherwise..
                continue

            # file key should go to result dict no matter what
            result[key] = filename
            logger.debug("Adding key %s: %s", key, filename)

            # filename resolved, check if exists and if given, hash
            file_path = self.raw_data.joinpath(filename)
            if file_path.is_file():
                if file_hash and file_hash != calculate_sha256(file_path):
                    logfunc(key, "Invalid hash for file %s!", filename)
            else:
                logfunc(
                    key,
                    "% referenced in %s but not found in archive!",
                    key,
                    METADATA_FILE,
                )

        files = []
        for f in list(required) + list(recommended):
            try:
                files.append(result[f])
            except KeyError:
                logfunc(
                    f,
                    "%: file %s expected but not found in %s file!",
                    obj_identifier,
                    f,
                    METADATA_FILE,
                )
                files.append(None)

        logger.debug("Returning files: %s", files)

        return files

    @create_objects(depth=1)
    def process_experiment(
        self, item_data: tuple[str, dict] | None = None, **kwargs
    ) -> ProcessedObject | None:
        """Extract data from yaml block for creating Experiment instance.

        Incoming data format (relevant bits):
        (
            protein_name: <str>,
            {
                'type': 'manual',
                'crystallographic_files': {
                    'xtal_pdb': {
                        'file': 'upload_1/crystallographic_files/5rgs/5rgs.pdb',
                        'sha256': sha <str>,
                    },
                    'xtal_mtz': {
                        'file': 'upload_1/crystallographic_files/5rgs/5rgs.mtz',
                        'sha256': sha <str>,
                    }
                },
                'status': 'new',
            }
        )

        This is enough to save full instance
        """
        logger.debug("incoming data: %s", item_data)
        experiment_name, data = item_data

        if "aligned_files" not in data.keys():
            return None

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="crystals",
            item_name=experiment_name,
        )

        (  # pylint: disable=unbalanced-tuple-unpacking
            pdb_info,
            mtz_info,
            cif_info,
        ) = self.validate_files(
            obj_identifier=experiment_name,
            file_struct=data["crystallographic_files"],
            required=("xtal_pdb", "xtal_mtz"),
            recommended=("ligand_cif",),
        )

        dtype = extract(key="type")

        if dtype == "manual":
            exp_type = 1
        elif dtype == "model_building":
            exp_type = 0
        else:
            exp_type = -1
            self.report.log(
                Level.FATAL, f"Unexpected 'type' '{dtype}' value for {experiment_name}"
            )

        dstatus = extract(key="status")

        if dstatus == "new":
            status = 0
        elif dstatus == "deprecated":
            status = 1
        elif dstatus == "superseded":
            status = 2
        else:
            status = -1
            self.report.log(
                Level.FATAL, f"Unexpected status '{dstatus}' for {experiment_name}"
            )

        # TODO: unhandled atm
        # version	int	old versions are kept	target loader
        version = 1

        fields = {
            "experiment_upload": self.experiment_upload,
            "code": experiment_name,
        }
        defaults = {
            "status": status,
            "version": version,
            "type": exp_type,
            "pdb_info": str(self._get_final_path(pdb_info)),
            "mtz_info": str(self._get_final_path(mtz_info)),
            "cif_info": str(self._get_final_path(cif_info)),
            # this doesn't seem to be present
            # pdb_sha256:
        }

        assigned_xtalform = extract(key="assigned_xtalform")

        index_fields = {
            "xtalform": assigned_xtalform,
        }

        return ProcessedObject(
            model_class=Experiment,
            fields=fields,
            key=experiment_name,
            defaults=defaults,
            index_data=index_fields,
        )

    @create_objects(depth=1)
    def process_compound(
        self, item_data: tuple[str, dict] | None = None, **kwargs
    ) -> ProcessedObject | None:
        """Extract data from yaml block for creating Compound instance.

        Incoming data format:
        xtal_pdb: {file: <file path>, sha256: <hash>}
        xtal_mtz: {file: <file path>, sha256: <hash>}
        ligand_cif: {file: <file path>, sha256: <hash>, smiles: <smiles>}
        panddas_event_files:
        - {file: <file path>, sha256: <hash>,
          model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}
        - {file: <file path>, sha256: <hash>,
          model: <int>, chain: <char[1]>, res: <int>, index: <int>, bdc: <float>}

        NB! After creation, many2many with project needs to be populated
        """
        logger.debug("incoming data: %s", item_data)
        protein_name, data = item_data
        if (
            "aligned_files" not in data.keys()
            or "crystallographic_files" not in data.keys()
        ):
            return None

        try:
            smiles = data["crystallographic_files"]["ligand_cif"]["smiles"]
        except KeyError as exc:
            # just setting the var to something
            smiles = (
                "crystallographic_files"
                if exc.args[0] == "ligand_cif"
                else "ligand_cif"
            )
            self.report.log(
                Level.FATAL,
                "{} missing from {} in '{}' experiment section".format(
                    exc, smiles, protein_name
                ),
            )

        fields = {
            "smiles": smiles,
        }

        return ProcessedObject(
            model_class=Compound,
            fields=fields,
            key=protein_name,
        )

    @create_objects(depth=1)
    def process_xtalform(
        self, item_data: tuple[str, dict] | None = None, **kwargs
    ) -> ProcessedObject | None:
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
        # weirdly, none of the fields is mandatory in Xtalform
        xtalform_name, data = item_data

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="xtalforms",
            item_name=xtalform_name,
        )

        fields = {
            "name": xtalform_name,
        }
        space_group = extract(key="xtalform_space_group")
        unit_cell_info = extract(key="xtalform_cell")

        defaults = {
            "space_group": space_group,
            "unit_cell_info": unit_cell_info,
        }

        return ProcessedObject(
            model_class=Xtalform,
            fields=fields,
            key=xtalform_name,
            defaults=defaults,
        )

    @create_objects(depth=1)
    def process_quat_assembly(
        self, item_data: tuple[str, dict] | None = None, **kwargs
    ) -> ProcessedObject | None:
        """Create QuatAssemblylform model instance from data.

        Incoming data format:
        <idx>:
            reference: <name>
            biomol: <biomol: str>
            chains: <chain info: str>

        No references to other models.
        """
        assembly_name, data = item_data

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="assemblies",
            item_name=assembly_name,
        )

        chains = extract(key="chains", level=Level.WARNING)

        fields = {
            "name": assembly_name,
            "chains": chains,
        }

        return ProcessedObject(
            model_class=QuatAssembly,
            fields=fields,
            key=assembly_name,
        )

    @create_objects(depth=3)
    def process_xtalform_quatassembly(
        self,
        xtalforms: dict[int | str, MetadataObject],
        quat_assemblies: dict[int | str, MetadataObject],
        item_data: tuple[str, str, int, dict] | None = None,
        **kwargs,
    ) -> ProcessedObject | None:
        """Create XtalformQuatAssembly model instance from data.

        Incoming data format:
        <idx>:
            assembly: <assembly id: int>
            chains: <str>

        """
        xtalform_id, _, _, data = item_data

        # hm.. doesn't reflect the fact that it's from a different
        # file.. and the message should perhaps be a bit different
        extract = functools.partial(
            self._extract,
            data=data,
            section_name="xtalforms",
            item_name=xtalform_id,
        )

        xtalform = xtalforms[xtalform_id].instance

        assembly_id = extract(key="assembly")

        # TODO: need to key check these as well..
        assembly = quat_assemblies[assembly_id].instance

        fields = {
            "assembly_id": assembly_id,
            "xtalform": xtalform,
            "quat_assembly": assembly,
            "chains": data["chains"],
        }

        return ProcessedObject(
            model_class=XtalformQuatAssembly,
            fields=fields,
            key=xtalform_id,
        )

    @create_objects(depth=1)
    def process_canon_site(
        self, item_data: tuple[str, dict] | None = None, **kwargs
    ) -> ProcessedObject | None:
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
        canon_site_id, data = item_data

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="canon_sites",
            item_name=canon_site_id,
        )

        residues = extract(key="residues", return_type=list)

        fields = {
            "name": canon_site_id,
            "residues": residues,
        }

        conf_sites_ids = extract(key="conformer_site_ids", return_type=list)
        ref_conf_site_id = extract(key="reference_conformer_site_id")

        index_data = {
            "ref_conf_site": ref_conf_site_id,
            "conformer_site_ids": conf_sites_ids,
            "reference_conformer_site_id": ref_conf_site_id,
        }

        return ProcessedObject(
            model_class=CanonSite,
            fields=fields,
            index_data=index_data,
            key=canon_site_id,
        )

    @create_objects(depth=1)
    def process_canon_site_conf(
        self,
        canon_sites: dict[str, ModelInstance],
        item_data: tuple[str, dict] | None = None,
        **kwargs,
    ) -> ProcessedObject | None:
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx: str>:
          reference_ligand_id: <lig_ref>
          residues: <array[char]>
          members: <array[char]>

        Unable to add references to:
        - SiteObservation (ref_site_observation)
        """
        conf_site_name, data = item_data
        canon_site = canon_sites[conf_site_name]

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="conformer_sites",
            item_name=conf_site_name,
            return_type=list,
        )

        residues = extract(key="residues")

        fields = {
            "name": conf_site_name,
            "residues": residues,
            "canon_site": canon_site,
        }

        # members = extract(key="members")
        ref_ligands = extract(key="reference_ligand_id")

        index_fields = {
            # "members": members,
            "reference_ligands": ref_ligands,
        }

        return ProcessedObject(
            model_class=CanonSiteConf,
            fields=fields,
            index_data=index_fields,
            key=conf_site_name,
        )

    @create_objects(depth=1)
    def process_xtalform_site(
        self,
        xtalforms: dict[int | str, MetadataObject],
        canon_sites: dict[str, ModelInstance],
        item_data: tuple[str, dict] | None = None,
        **kwargs,
    ) -> ProcessedObject | None:
        """Create Xtalform model instance from data.

        Incoming data format:
        <idx>:
          xtalform_id: <str>
          canonical_site_id: <str>
          crystallographic_chain: A
          members: <array[str]>

        Saves references to all other tables (Xtalform and CanonSite).
        """
        xtalform_site_name, data = item_data

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="xtalform_sites",
            item_name=xtalform_site_name,
        )

        xtalform_id = extract(key="xtalform_id")

        canon_site_id = extract(key="canonical_site_id")

        xtalform = xtalforms[xtalform_id].instance
        canon_site = canon_sites[canon_site_id]

        lig_chain = extract(key="crystallographic_chain")
        residues = extract(key="members", return_type=list)

        fields = {
            "xtalform_site_id": xtalform_site_name,
            "xtalform": xtalform,
            "canon_site": canon_site,
        }

        defaults = {
            "lig_chain": lig_chain,
            "residues": residues,
        }

        return ProcessedObject(
            model_class=XtalformSite,
            fields=fields,
            defaults=defaults,
            key=xtalform_site_name,
        )

    @create_objects(depth=5)
    def process_site_observation(
        self,
        experiments: dict[int | str, MetadataObject],
        compounds: dict[int | str, MetadataObject],
        xtalform_sites: dict[str, ModelInstance],
        canon_site_confs: dict[int | str, MetadataObject],
        item_data: tuple[str, str, str, int | str, int | str, dict] | None = None,
        # chain: str,
        # ligand: str,
        # idx: int | str,
        # data: dict,
        **kwargs,
    ) -> ProcessedObject | None:
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
        try:
            experiment_id, _, chain, ligand, idx, data = item_data
        except ValueError:
            # wrong data item
            return None

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="crystals",
            item_name=experiment_id,
            level=Level.INFO,
        )

        experiment = experiments[experiment_id].instance

        code = f"{experiment.code}_{chain}_{str(ligand)}_{str(idx)}"
        key = f"{experiment.code}/{chain}/{str(ligand)}"

        compound = compounds[experiment_id].instance
        canon_site_conf = canon_site_confs[idx].instance
        xtalform_site = xtalform_sites[key]

        (  # pylint: disable=unbalanced-tuple-unpacking
            bound_file,
            apo_solv_file,
            apo_desolv_file,
            apo_file,
            xmap_2fofc_file,
            xmap_fofc_file,
            event_file,
            artefacts_file,
            ligand_mol,
        ) = self.validate_files(
            obj_identifier=experiment_id,
            file_struct=data,
            required=(
                "structure",
                "pdb_apo_solv",
                "pdb_apo_desolv",
                "pdb_apo",
                "2Fo-Fc_map",
                "Fo-Fc_map",
                "event_map",
            ),
            recommended=(
                "artefacts",
                "ligand_mol",
            ),
        )

        mol_data = None
        try:
            with open(
                self.raw_data.joinpath(ligand_mol),
                "r",
                encoding="utf-8",
            ) as f:
                mol_data = f.read()
        except TypeError:
            # this site observation doesn't have a ligand. perfectly
            # legitimate case
            pass

        smiles = extract(key="ligand_smiles")

        fields = {
            # Code for this protein (e.g. Mpro_Nterm-x0029_A_501_0)
            "code": code,
            "experiment": experiment,
            "cmpd": compound,
            "xtalform_site": xtalform_site,
            "canon_site_conf": canon_site_conf,
            "smiles": smiles,
            "seq_id": ligand,
            "chain_id": chain,
        }

        defaults = {
            "bound_file": str(self._get_final_path(bound_file)),
            "apo_solv_file": str(self._get_final_path(apo_solv_file)),
            "apo_desolv_file": str(self._get_final_path(apo_desolv_file)),
            "apo_file": str(self._get_final_path(apo_file)),
            "xmap_2fofc_file": str(self._get_final_path(xmap_2fofc_file)),
            "xmap_fofc_file": str(self._get_final_path(xmap_fofc_file)),
            "event_file": str(self._get_final_path(event_file)),
            "artefacts_file": str(self._get_final_path(artefacts_file)),
            "pdb_header_file": "currently missing",
            "ligand_mol_file": mol_data,
        }

        return ProcessedObject(
            model_class=SiteObservation,
            fields=fields,
            defaults=defaults,
            key=key,
        )

    def process_metadata(
        self,
        upload_root: Path,
    ):
        """Extract model instances from metadata file and save them to db."""
        # TODO: this method is quite long and should perhaps be broken
        # apart. then again, it just keeps calling other methods to
        # create model instances, so logically it's just doing the
        # same thing.

        # update_task(task, "PROCESSING", "Processing metadata")
        logger.info("%sProcessing %s", self.report.task_id, upload_root)

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

        (  # pylint: disable=unbalanced-tuple-unpacking
            assemblies,
            xtalform_assemblies,
        ) = self._get_yaml_blocks(
            yaml_data=self._load_yaml(Path(upload_root).joinpath(XTALFORMS_FILE)),
            blocks=("assemblies", "xtalforms"),
        )

        meta = self._load_yaml(Path(upload_root).joinpath(METADATA_FILE))

        # collect top level info
        self.version_number = meta["version_number"]
        self.version_dir = meta["version_dir"]
        self.previous_version_dirs = meta["previous_version_dirs"]

        (  # pylint: disable=unbalanced-tuple-unpacking
            crystals,
            xtalforms,
            canon_sites,
            conformer_sites,
            xtalform_sites,
        ) = self._get_yaml_blocks(
            yaml_data=meta,
            blocks=(
                "crystals",
                "xtalforms",
                "canon_sites",
                "conformer_sites",
                "xtalform_sites",
            ),
        )

        experiment_objects = self.process_experiment(yaml_data=crystals)
        compound_objects = self.process_compound(yaml_data=crystals)

        # save components manytomany to experiment
        # TODO: is it 1:1 relationship? looking at the meta_align it
        # seems to be, but why the m2m then?
        for (
            comp_code,
            comp_meta,
        ) in compound_objects.items():  # pylint: disable=no-member
            experiment = experiment_objects[comp_code].instance
            experiment.compounds.add(comp_meta.instance)

        xtalform_objects = self.process_xtalform(yaml_data=xtalforms)

        # add xtalform fk to experiment
        for _, obj in experiment_objects.items():  # pylint: disable=no-member
            try:
                obj.instance.xtalform = xtalform_objects[
                    obj.index_data["xtalform"]
                ].instance
                obj.instance.save()
            except KeyError:
                # TODO: message may need tweaking
                msg = f"xtalform {obj.instance.code} undefined for {obj}"
                logger.warning(msg)

        quat_assembly_objects = self.process_quat_assembly(yaml_data=assemblies)

        _ = self.process_xtalform_quatassembly(
            yaml_data=xtalform_assemblies,
            xtalforms=xtalform_objects,
            quat_assemblies=quat_assembly_objects,
        )

        canon_site_objects = self.process_canon_site(yaml_data=canon_sites)
        # NB! missing fk's:
        # - ref_conf_site
        # - quat_assembly

        # reindex canon sites by canon_sites_conf_sites
        # NB! this is also used below for ref_conf_site in canon_site
        canon_sites_by_conf_sites = {
            conf: obj.instance
            for obj in canon_site_objects.values()  # pylint: disable=no-member
            for conf in obj.index_data["conformer_site_ids"]
        }

        canon_site_conf_objects = self.process_canon_site_conf(
            yaml_data=conformer_sites, canon_sites=canon_sites_by_conf_sites
        )
        # NB! missing fk's:
        # - site_observation

        xtalform_sites_objects = self.process_xtalform_site(
            yaml_data=xtalform_sites,
            canon_sites=canon_sites_by_conf_sites,
            xtalforms=xtalform_objects,
        )

        # now can update CanonSite with ref_conf_site
        for val in canon_site_objects.values():  # pylint: disable=no-member
            val.instance.ref_conf_site = canon_site_conf_objects[
                val.index_data["reference_conformer_site_id"]
            ].instance
            val.instance.save()

        # canon site instances are now complete
        # still missing fk to site_observation in canon_site_conf

        # reindex xtalform site to grab for site observation I don't
        # need this anywhere else, why won't i just give the correct
        # key for xtal sites objects?
        xtalform_site_by_tag = {}
        for val in xtalform_sites_objects.values():  # pylint: disable=no-member
            for k in val.instance.residues:
                xtalform_site_by_tag[k] = val.instance

        site_observation_objects = self.process_site_observation(
            yaml_data=crystals,
            experiments=experiment_objects,
            compounds=compound_objects,
            xtalform_sites=xtalform_site_by_tag,
            canon_site_confs=canon_site_conf_objects,
        )

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
        for val in canon_site_conf_objects.values():  # pylint: disable=no-member
            val.instance.ref_site_observation = site_observation_objects[
                val.index_data["reference_ligands"]
            ].instance
            val.instance.save()

    def process_bundle(self):
        """Resolves subdirs in uploaded data bundle.

        If called from task, takes task as a parameter for status updates.
        """

        # by now I should have archive unpacked, get target name from
        # config.yaml
        up_iter = self.raw_data.glob("upload_*")
        try:
            upload_dir = next(up_iter)
        except StopIteration as exc:
            msg = "Upload directory missing from uploaded file!"
            self.report.log(Level.FATAL, msg)
            raise StopIteration(msg) from exc

        try:
            upload_dir = next(up_iter)
            self.report.log(Level.WARNING, "Multiple upload directories in archive!")
        except StopIteration as exc:
            # just a warning, ignoring the second one
            pass

        # now that target name is not included in path, I don't need
        # it here, need it just before creating target object. Also,
        # there's probably no need to throw a fatal here, I can
        # reasonably well deduce it from meta (I think)
        config_it = upload_dir.glob(CONFIG_FILE)
        try:
            config_file = next(config_it)
        except StopIteration as exc:
            raise StopIteration(f"config file missing from {str(upload_dir)}") from exc

        config = self._load_yaml(config_file)
        logger.debug("config: %s", config)

        try:
            self.target_name = config["target_name"]
        except KeyError as exc:
            raise KeyError("target_name missing in config file!") from exc

        try:
            self.process_metadata(
                upload_root=upload_dir,
            )
        except FileNotFoundError as exc:
            # result.append(exc.args[0])
            raise FileNotFoundError(exc.args[0]) from exc
        except IntegrityError as exc:
            # result.append(exc.args[0])
            raise IntegrityError(exc.args[0]) from exc
        except ValueError as exc:
            # result.append(exc.args[0])
            raise ValueError(exc.args[0]) from exc
        except FileExistsError as exc:
            # this target has been uploaded at- least once and
            # this upload already exists. skip
            # result.append(exc.args[0])
            raise FileExistsError(exc.args[0]) from exc

    def _load_yaml(self, yaml_file: Path) -> dict:
        try:
            with open(yaml_file, "r", encoding="utf-8") as file:
                contents = yaml.safe_load(file)
        except FileNotFoundError as exc:
            msg = f"{yaml_file.stem} file not found in data archive!"
            # logger.error("%s%s", self.task_id, msg)
            raise FileNotFoundError(msg) from exc

        return contents

    # TODOL error handling. what's the correct response when
    # something's missing? push through and compile report?
    def _get_yaml_blocks(self, yaml_data: dict, blocks: Iterable) -> list[dict]:
        validation_errors = []
        error_text = "'{}' section missing in input file"
        result = []
        for block in blocks:
            try:
                result.append(yaml_data[block])
            except KeyError:
                msg = error_text.format(block)
                validation_errors.append(msg)
                # update_task(self.task, "ERROR", msg)
                # logger.error("%s%s", self.task_id, msg)

        # wonder if it's worth throwing a custom exception here.. easier
        # to scan the logs. then again, if errors are passed to user,
        # inspecting logs isn't really necessary?
        if validation_errors:
            raise KeyError(validation_errors)

        return result

    def _extract(
        self,
        data: dict,
        key: str | int,
        section_name: str,
        item_name: str,
        level: Level = Level.FATAL,
        return_type: type = str,
    ) -> Any:
        try:
            result = data[key]
        except KeyError as exc:
            if level == Level.INFO:
                result = ""
            else:
                result = "missing"
            if return_type == list:
                result = [result]

            self.report.log(
                level,
                "{} missing from {}: {} section".format(
                    exc,
                    section_name,
                    item_name,
                ),
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

    def _get_final_path(self, path: str | None) -> Path | None:
        """Update relative path to final storage path

        NB! this returns a relative path that can be used in queries
        not absoulte one. This is used to populate location fields in
        database tables.
        """
        try:
            return self.final_path.joinpath(path)
        except TypeError:
            # received invalid path
            return None


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
                logger.info("%s%s", target_loader.report.task_id, msg)
                # update_task(task, "PROCESSING", msg)
                archive.extractall(target_loader.raw_data)
                msg = f"Data extraction complete: {data_bundle}"
                logger.info("%s%s", target_loader.report.task_id, msg)
        except FileNotFoundError as exc:
            msg = f"{data_bundle} file does not exist!"
            logger.exception("%s%s", target_loader.report.task_id, msg)
            target_loader.experiment_upload.message = exc.args[0]
            raise FileNotFoundError(msg) from exc

        try:
            with transaction.atomic():
                target_loader.process_bundle()
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

        set_directory_permissions(target_loader.abs_final_path, 0o755)

        target_loader.report.final()
        target_loader.experiment_upload.message = target_loader.report.json()

        # logger.debug("%s", upload_report)

        target_loader.experiment_upload.save()
