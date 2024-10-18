import contextlib
import functools
import hashlib
import logging
import math
import os
import shutil
import tarfile
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Iterable, List, Optional, Tuple, TypeVar

import pandas as pd
import yaml
from celery import Task
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.postgres.aggregates import ArrayAgg
from django.core.exceptions import MultipleObjectsReturned
from django.db import IntegrityError, transaction
from django.db.models import Count, F, Model
from django.db.models.base import ModelBase
from django.utils import timezone

from api.utils import deployment_mode_is_production
from fragalysis.settings import TARGET_LOADER_MEDIA_DIRECTORY
from scoring.models import SiteObservationGroup
from viewer.models import (
    CanonSite,
    CanonSiteConf,
    Compound,
    CompoundIdentifier,
    CompoundIdentifierType,
    Experiment,
    ExperimentUpload,
    Pose,
    Project,
    QuatAssembly,
    SiteObservation,
    SiteObservationTag,
    TagCategory,
    Target,
    Xtalform,
    XtalformQuatAssembly,
    XtalformSite,
)
from viewer.utils import alphanumerator, clean_object_id, sanitize_directory_name

logger = logging.getLogger(__name__)

# data that goes to tables are in the following files
# assemblies and xtalforms
XTALFORMS_FILE = "assemblies.yaml"

# target name, nothing else
CONFIG_FILE = "config*.yaml"

# everything else
METADATA_FILE = "meta_aligner.yaml"

# transformation matrices
TRANS_NEIGHBOURHOOD = "neighbourhood_transforms.yaml"
TRANS_CONF_SITE = "conformer_site_transforms.yaml"
TRANS_REF_STRUCT = "reference_structure_transforms.yaml"

CUSTOM_IDENTIFIER_FILE = "compounds_manual.csv"
CUSTOM_IDENTIFIER_COLUMNS = ("xtal", "ligand_name", "compound_code")


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
    CANCELED = "CANCELED"


@dataclass
class MetadataObject:
    """Data structure to store freshly created model instances.

    data blocks from from meta_aligner.yaml are processed into
    dictionaries: { some_id: MetadataObjects, ...}

    Reason being, quite often I need to refer to these by some
    alternative ID. With the dataclass, I'm able to create temporary
    dicts with key that are needed.
    """

    instance: Model
    key: str
    versioned_key: str
    index_data: dict = field(default_factory=dict)
    new: bool = False


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
    key: str | tuple[str, str]
    defaults: dict = field(default_factory=dict)
    index_data: dict = field(default_factory=dict)
    versioned_key: Optional[str | tuple[str, str]] = ""


@dataclass
class UploadReportEntry:
    message: str
    level: int | None = None

    def __str__(self):
        if self.level is None:
            return self.message
        return f"{logging.getLevelName(self.level)}: {self.message}"


@dataclass
class UploadReport:
    task: Task | None
    proposal_ref: str
    stack: list[UploadReportEntry] = field(default_factory=list)
    upload_state: UploadState = UploadState.PROCESSING
    failed: bool = False

    def __post_init__(self) -> None:
        self.task_id = f"task {self.task.request.id}: " if self.task else ""

    def log(self, level: int, message: str) -> None:
        msg = f"{self.task_id}{message}"
        if level == logging.ERROR:
            self.failed = True
            self.upload_state = UploadState.REPORTING
        logger.log(level, msg)
        self.stack.append(UploadReportEntry(level=level, message=message))
        self._update_task(self.json())

    def final(self, message, success=True):
        self.upload_state = UploadState.SUCCESS

        # This is (expected to be) the last message for the upload.
        # Add the user-supplied message and then add a string indicating success or failure.
        self.stack.append(UploadReportEntry(message=message))
        status_line = 'SUCCESS' if success else 'FAILED'
        self.stack.append(UploadReportEntry(message=status_line))

        self._update_task(self.json())

    def json(self):
        return [str(k) for k in self.stack]

    def _update_task(self, message: str | list) -> None:
        if not self.task:
            return
        with contextlib.suppress(AttributeError):
            self.task.update_state(
                state=self.upload_state,
                meta={
                    "proposal_ref": self.proposal_ref,
                    "description": message,
                },
            )


def _validate_bundle_against_mode(config_yaml: Dict[str, Any]) -> Optional[str]:
    """Inspects the meta to ensure it is supported by the MODE this stack is in.
    Mode is (typically) one of DEVELOPER or PRODUCTION.
    """
    assert config_yaml
    if not deployment_mode_is_production():
        # We're not in production mode - no bundle checks
        return None

    # PRODUCTION mode (strict)

    # Initial concern - the loader's git information.
    # It must not be 'dirty' and must have a valid 'tag'.
    xca_git_info_key = "xca_git_info"
    base_error_msg = "Stack is in PRODUCTION mode - and"
    try:
        xca_git_info = config_yaml[xca_git_info_key]
    except KeyError:
        return f"{base_error_msg} '{xca_git_info_key}' is a required configuration property"

    logger.info("%s: %s", xca_git_info_key, xca_git_info)

    if "dirty" not in xca_git_info:
        return f"{base_error_msg} '{xca_git_info_key}' has no 'dirty' property"
    if xca_git_info["dirty"]:
        return f"{base_error_msg} '{xca_git_info_key}->dirty' must be False"

    if "tag" not in xca_git_info:
        return f"{base_error_msg} '{xca_git_info_key}' has no 'tag' property"
    xca_version_tag: str = str(xca_git_info["tag"])
    tag_parts: List[str] = xca_version_tag.split(".")
    tag_valid: bool = True
    if len(tag_parts) in {2, 3}:
        for tag_part in tag_parts:
            if not tag_part.isdigit():
                tag_valid = False
                break
    else:
        tag_valid = False
    if not tag_valid:
        return f"{base_error_msg} '{xca_git_info_key}->tag' must be 'N.N[.N]'. Got '{xca_version_tag}'"

    # OK if we get here
    return None


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


def set_directory_permissions(path, permissions) -> None:
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
def calculate_sha256(filepath) -> str:
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read the file in chunks of 4096 bytes
        for chunk in iter(lambda: f.read(4096), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def strip_version(s: str, separator: str = "/") -> Tuple[str, int]:
    # format something like XX01ZVNS2B-x0673/B/501/1
    # remove tailing '<separator>1'
    return s[0 : s.rfind(separator)], int(s[s.rfind(separator) + 1 :])


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
    ) -> dict[int | str | tuple[str, str], MetadataObject]:
        # logger.debug("+ wrapper_service_query")
        # logger.debug("args passed: %s", args)
        # logger.debug("kwargs passed: %s", kwargs)

        flattened_data = flatten_dict(yaml_data, depth=depth)
        result = {}
        created, existing, failed, updated = 0, 0, 0, 0

        for item in flattened_data:
            logger.debug("flattened data item: %s", item)
            instance_data = func(
                self, *args, item_data=item, validate_files=False, **kwargs
            )
            logger.debug("Instance data returned: %s", instance_data)
            obj = None
            new = False
            if not instance_data:
                continue

            try:
                if instance_data.fields:
                    try:
                        obj = instance_data.model_class.filter_manager.by_target(
                            self.target
                        ).get(**instance_data.fields)
                        logger.debug("Object exists: %s", instance_data.fields)
                        new = False
                    except instance_data.model_class.DoesNotExist:
                        # revalidate files
                        logger.debug("Object doesn't exist: %s", instance_data)
                        instance_data = func(self, *args, item_data=item, **kwargs)
                        obj = instance_data.model_class(
                            **instance_data.fields,
                            **instance_data.defaults,
                        )
                        obj.save()
                        new = True
                    except MultipleObjectsReturned:
                        msg = "{}.get_or_create in {} returned multiple objects for {}".format(
                            instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                            instance_data.key,
                            instance_data.fields,
                        )
                        self.report.log(logging.ERROR, msg)
                        failed = failed + 1

                else:
                    # no unique field requirements, just create new object
                    obj = instance_data.model_class(
                        **instance_data.defaults,
                    )
                    obj.save()
                    new = True
                logger.debug(
                    "%s object %s created",
                    instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                    obj,
                )
            except IntegrityError:
                msg = "{} object {} failed to save".format(
                    instance_data.model_class._meta.object_name,  # pylint: disable=protected-access
                    instance_data.key,
                )
                self.report.log(logging.ERROR, msg)
                failed = failed + 1

            if obj:
                # update any additional fields
                instance_qs = instance_data.model_class.objects.filter(pk=obj.pk)
                instance_qs.update(**instance_data.defaults)
                obj.refresh_from_db()
            else:
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

            if new:
                created = created + 1
                # check if old versions exist and mark them as superseded
                if "version" in instance_data.fields.keys():
                    del instance_data.fields["version"]
                    superseded = instance_data.model_class.objects.filter(
                        **instance_data.fields,
                    ).exclude(
                        pk=obj.pk,
                    )
                    updated += superseded.update(superseded=True)

            else:
                existing = existing + 1

            m = MetadataObject(
                instance=obj,
                key=instance_data.key,
                versioned_key=instance_data.versioned_key,
                index_data=instance_data.index_data,
                new=new,
            )
            # index data here probs
            result[instance_data.versioned_key] = m

        if result:
            msg = "{} {} objects processed, {} created, {} fetched from database".format(
                created + existing + failed,
                next(  # pylint: disable=protected-access
                    iter(result.values())
                ).instance._meta.model._meta.object_name,  # pylint: disable=protected-access
                created,
                existing,
            )  # pylint: disable=protected-access
            self.report.log(logging.INFO, msg)
        else:
            # cannot continue when one object type is missing, abort
            msg = f"No objects returned by {func.__name__}"
            self.report.log(logging.ERROR, msg)

        # refresh all objects to make sure they're up to date.
        # this is specifically because of the superseded flag above -
        # I'm setting this in separate queryset, the db rows are
        # updated, but the changes are not being propagated to the
        # objects in result dict. Well aware that this isn't efficient
        # but I don't have access to parent's versioned key here, (and
        # even if I did, there's no guarantee that they would have
        # already been processd), so that's why updating every single
        # object.
        if updated > 0:
            for k in result.values():
                k.instance.refresh_from_db()

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
        self.version_number = 1
        self.version_dir = None
        self.previous_version_dirs = None
        self.user_id = user_id

        self.report = UploadReport(task=task, proposal_ref=self.proposal_ref)

        self.raw_data.mkdir()

        # create exp upload object
        # NB! this is not saved here in case upload fails
        self.experiment_upload = ExperimentUpload(
            commit_datetime=timezone.now(),
            file=self.data_bundle,
        )

        # work out where the data finally lands
        path = Path(TARGET_LOADER_MEDIA_DIRECTORY)

        # give each upload a unique directory
        # update: resolving issue 1311 introduced a bug, where
        # subsequent uploads overwrote file paths and files appeared
        # to be missing. changing the directory structure so this
        # wouldn't be an issue, the new structure is
        # target_loader_data/target_title/upload_(n)/...
        if task:
            self.experiment_upload.task_id = task.request.id

        # figure out absolute and relative paths to final
        # location. relative path is added to db field, this will be
        # used in url requests to retrieve the file. absolute path is
        # for moving the file to the final location
        self._final_path = path
        self._abs_final_path = Path(settings.MEDIA_ROOT).joinpath(path)
        # but don't create now, this comes later

        # to be used in logging messages, if no task, means invoked
        # directly, likely from management command
        # self.task_id = f"task {task.request.id}: " if task else ""

        # these will be filled later
        self.target_name = None
        self._target_root = None
        self.target = None
        self.project = None

        # Initial (reassuring message)
        bundle_filename = os.path.basename(self.bundle_path)
        self.report.log(
            logging.INFO,
            f"Created TargetLoader for '{bundle_filename}' proposal_ref='{proposal_ref}'",
        )

    @property
    def final_path(self) -> Path:
        return self._final_path

    @property
    def abs_final_path(self) -> Path:
        return self._abs_final_path

    def validate_map_files(
        self,
        key: str,
        obj_identifier: str,
        file_struct: list,
        validate_files: bool = True,
    ) -> tuple[list[str], list[str]]:
        """Validate list of panddas event files.

        Special case of file validation, too complex to squeeze into
        the main validation method (mainly because of typing).
        """

        def logfunc(_, message):
            self.report.log(logging.WARNING, message)

        paths = []
        source_files = []
        for item in file_struct:
            fname, file_hash = self._check_file(item, obj_identifier, key, logfunc)
            source_file = item.get("source_file", None)
            if not fname:
                continue

            if validate_files:
                self._check_file_hash(obj_identifier, key, fname, file_hash, logfunc)
            paths.append(fname)
            source_files.append(source_file)

        return paths, source_files

    def validate_files(
        self,
        obj_identifier: str,
        file_struct: dict,
        required: Iterable[str] = (),
        recommended: Iterable[str] = (),
        validate_files: bool = True,
    ) -> list[tuple[str | None, str | None]]:
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

        def logfunc(key, message):
            if key in required:
                self.report.log(logging.ERROR, message)
            else:
                self.report.log(logging.WARNING, message)

        result = {}
        for key, value in file_struct.items():
            if key not in required and key not in recommended:
                # schema isn't looking for this file, ignore
                continue

            filename, file_hash, source_file = None, None, None

            # sort out the filename
            if isinstance(value, dict):
                filename, file_hash = self._check_file(
                    value, obj_identifier, key, logfunc
                )
                if not filename:
                    continue

                source_file = value.get("source_file", None)

                if validate_files:
                    self._check_file_hash(
                        obj_identifier, key, filename, file_hash, logfunc
                    )

            elif isinstance(value, str):
                filename = value
                if validate_files:
                    self._check_file_hash(
                        obj_identifier, key, filename, file_hash, logfunc
                    )

            else:
                # probably panddas files here
                continue

            # file key should go to result dict no matter what
            result[key] = (filename, source_file)
            logger.debug("Adding key %s: %s", key, filename)

        files = []
        for f in list(required) + list(recommended):
            try:
                files.append(result[f])
            except KeyError:
                logfunc(
                    f,
                    "{}: file {} expected but not found in {} file".format(
                        obj_identifier,
                        f,
                        METADATA_FILE,
                    ),
                )
                files.append((None, None))  # type: ignore [arg-type]

        logger.debug("Returning files: %s", files)

        # memo to self: added type ignore directives to return line
        # below and append line above because after small refactoring,
        # mypy all of the sudden started throwing errors on both of
        # these. the core of it's grievance is that it expects the
        # return type to be list[str]. no idea why, function signature
        # clearly defines it as list[str | None]

        return files  # type: ignore [return-value]

    def _check_file(
        self,
        value: dict,
        obj_identifier: str,
        key: str,
        logfunc: Callable,
    ) -> Tuple[str | None, str | None]:
        file_hash = value.get("sha256")
        try:
            filename = value["file"]
        except KeyError:
            # this is rather unexpected, haven't seen it yet
            filename = None
            logfunc(key, f"{obj_identifier}: malformed dict, key 'file' missing")
        return filename, file_hash

    def _check_file_hash(
        self,
        obj_identifier: str,
        key: str,
        filename: str,
        file_hash: str | None,
        logfunc: Callable,
    ) -> None:
        file_path = self.raw_data.joinpath(filename)
        if file_path.is_file():
            if file_hash and file_hash != calculate_sha256(file_path):
                logfunc(key, f"Invalid hash for file {filename}")
        else:
            logfunc(
                key,
                f"{key} referenced in {METADATA_FILE}: {obj_identifier} but not found in archive",
            )

    def _enumerate_objects(self, objects: dict, attr: str) -> None:
        # don't overwrite values already in database, get the current
        # max value and continue from there
        max_existing = 0
        for val in objects.values():  # pylint: disable=no-member
            value = getattr(val.instance, attr, 0)
            if value:
                max_existing = max(value, max_existing)

        if not max_existing:
            max_existing = 0

        for val in objects.values():  # pylint: disable=no-member
            value = getattr(val.instance, attr)
            if not value:
                max_existing = max_existing + 1
                setattr(val.instance, attr, max_existing)
                val.instance.save()

    @create_objects(depth=1)
    def process_experiment(
        self,
        item_data: tuple[str, dict] | None = None,
        prefix_tooltips: dict[str, str] | None = None,
        validate_files: bool = True,
        **kwargs,
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
                    },
                    'panddas_event_files': {
                        'file': <path>.ccp4,
                        'sha256': sha <str>,
                        'model': '1', chain: B, res: 203, index: 1, bdc: 0.23
                    },
                'status': 'new',
                },

            }
        )

        This is enough to save full instance
        """
        del kwargs
        assert item_data
        logger.debug("incoming data: %s", item_data)
        experiment_name, data = item_data

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="crystals",
            item_name=experiment_name,
        )

        (  # pylint: disable=unbalanced-tuple-unpacking
            pdb_info_t,
            mtz_info_t,
            cif_info_t,
        ) = self.validate_files(
            obj_identifier=experiment_name,
            file_struct=data["crystallographic_files"],
            recommended=(
                "xtal_pdb",
                "xtal_mtz",
                "ligand_cif",
            ),
            validate_files=validate_files,
        )

        pdb_info, pdb_info_source_file = pdb_info_t
        mtz_info, mtz_info_source_file = mtz_info_t
        cif_info, cif_info_source_file = cif_info_t

        try:
            event_files = data["crystallographic_files"]["ligand_binding_events"]
        except KeyError:
            event_files = []

        map_info_files, map_info_source_files = self.validate_map_files(
            key="ligand_binding_events",
            obj_identifier=experiment_name,
            file_struct=event_files,
            validate_files=validate_files,
        )

        dtype = extract(key="type")

        if dtype == "manual":
            exp_type = 1
        elif dtype == "model_building":
            exp_type = 0
        else:
            exp_type = -1
            self.report.log(
                logging.ERROR,
                f"Unexpected 'type' '{dtype}' value for {experiment_name}",
            )

        dstatus = extract(key="status")

        status_codes = {
            "new": 0,
            "deprecated": 1,
            "superseded": 2,
            "unchanged": 3,
        }

        try:
            status = status_codes[dstatus]
        except KeyError:
            status = -1
            self.report.log(
                logging.ERROR, f"Unexpected status '{dstatus}' for {experiment_name}"
            )

        # if empty or key missing entirely, ensure code_prefix returns empty
        code_prefix = extract(key="code_prefix", level=logging.INFO)
        # ignoring type because tooltip dict can legitimately be empty
        # and in such case, assert statement fails. need to remove it
        # and use the ignore
        prefix_tooltip = prefix_tooltips.get(code_prefix, "")  # type: ignore[union-attr]

        fields = {
            "code": experiment_name,
        }

        map_info_paths = []
        if map_info_files:
            map_info_paths = [str(self._get_final_path(k)) for k in map_info_files]

        defaults = {
            "experiment_upload": self.experiment_upload,
            "status": status,
            "type": exp_type,
            "pdb_info": str(self._get_final_path(pdb_info)),
            "mtz_info": str(self._get_final_path(mtz_info)),
            "cif_info": str(self._get_final_path(cif_info)),
            "pdb_info_source_file": pdb_info_source_file,
            "mtz_info_source_file": mtz_info_source_file,
            "cif_info_source_file": cif_info_source_file,
            "map_info": map_info_paths,
            "map_info_source_files": map_info_source_files,
            "prefix_tooltip": prefix_tooltip,
            "code_prefix": code_prefix,
            # this doesn't seem to be present
            # pdb_sha256:
        }

        assigned_xtalform = extract(key="assigned_xtalform")

        index_fields = {
            "xtalform": assigned_xtalform,
            "code_prefix": code_prefix,
        }

        return ProcessedObject(
            model_class=Experiment,
            fields=fields,
            key=experiment_name,
            versioned_key=experiment_name,
            defaults=defaults,
            index_data=index_fields,
        )

    @create_objects(depth=5)
    def process_compound(
        self,
        experiments: dict[int | str, MetadataObject],
        item_data: tuple[str, str, str, str, str, dict] | None = None,
        **kwargs,
    ) -> ProcessedObject | None:
        """Extract data from yaml block for creating Compound instance.

        Incoming item_data format:
        experiment_name <str>,
        "crystallographic_files" <str>,
        "ligand_cif" <str>,
        "ligands" <str>,
        ligand_key <str>,
        {
            "smiles": smiles <str>
        }

        NB! After creation, many2many with project needs to be populated
        """
        del kwargs
        assert item_data
        logger.debug("incoming data: %s", item_data)

        # remove non-compound objects
        try:
            experiment_name, _, _, _, ligand_key, data = item_data
        except ValueError:
            # wrong data item
            return None

        # more validation
        if (
            item_data[1] != "crystallographic_files"
            or item_data[2] != "ligand_cif"
            or item_data[3] != "ligands"
            or not experiments[experiment_name].new  # remove already saved objects
        ):
            return None

        smiles = data.get("smiles", None)
        compound_code = data.get("compound_code", None)

        if smiles is None and compound_code is None:
            # gotta have at least something
            return None

        defaults = {
            "smiles": smiles,
            "compound_code": compound_code,
            "ligand_name": ligand_key,
        }

        return ProcessedObject(
            model_class=Compound,
            fields={},
            defaults=defaults,
            key=(experiment_name, ligand_key),
            versioned_key=(experiment_name, ligand_key),
        )

    @create_objects(depth=1)
    def process_xtalform(
        self,
        item_data: tuple[str, dict] | None = None,
        **kwargs,
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
        del kwargs
        assert item_data
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
            versioned_key=xtalform_name,
            defaults=defaults,
        )

    @create_objects(depth=1)
    def process_quat_assembly(
        self,
        item_data: tuple[str, dict] | None = None,
        **kwargs,
    ) -> ProcessedObject | None:
        """Create QuatAssemblylform model instance from data.

        Incoming data format:
        <idx>:
            reference: <name>
            biomol: <biomol: str>
            chains: <chain info: str>

        No references to other models.
        """
        del kwargs
        assert item_data
        assembly_name, data = item_data

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="assemblies",
            item_name=assembly_name,
        )

        chains = extract(key="chains", level=logging.WARNING)

        fields = {
            "name": assembly_name,
            "chains": chains,
        }

        return ProcessedObject(
            model_class=QuatAssembly,
            fields=fields,
            key=assembly_name,
            versioned_key=assembly_name,
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
        del kwargs
        assert item_data
        xtalform_id, _, assembly_id, data = item_data

        # hm.. doesn't reflect the fact that it's from a different
        # file.. and the message should perhaps be a bit different
        extract = functools.partial(
            self._extract,
            data=data,
            section_name="xtalforms",
            item_name=xtalform_id,
        )

        xtalform = xtalforms[xtalform_id].instance

        quat_assembly_id = extract(key="assembly")

        # TODO: need to key check these as well..
        assembly = quat_assemblies[quat_assembly_id].instance

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
            versioned_key=xtalform_id,
        )

    @create_objects(depth=1)
    def process_canon_site(
        self,
        item_data: tuple[str, dict] | None = None,
        **kwargs,
    ) -> ProcessedObject | None:
        """Create CanonSite model instance from data.

        Incoming data format:
        <id: str>:
            centroid_res: <str>
            conformer_site_ids: <array[str]>
            global_reference_dtag: <str>
            reference_conformer_site_id: <str>
            residues: <array[str]>

        Unable to add references to:
        - CanonSiteConf (ref_conf_site)

        """
        del kwargs
        assert item_data
        v_canon_site_id, data = item_data

        canon_site_id, version = strip_version(v_canon_site_id, separator="+")

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="canon_sites",
            item_name=canon_site_id,
        )

        residues = extract(key="residues", return_type=list)
        centroid_res = extract(key="centroid_res")
        conf_sites_ids = extract(key="conformer_site_ids", return_type=list)
        ref_conf_site_id = extract(key="reference_conformer_site_id")

        centroid_res = f"{centroid_res}_v{version}"

        fields = {
            "name": canon_site_id,
            "version": version,
        }

        defaults = {
            "residues": residues,
            "centroid_res": centroid_res,
        }

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
            versioned_key=v_canon_site_id,
            defaults=defaults,
        )

    @create_objects(depth=1)
    def process_canon_site_conf(
        self,
        canon_sites: dict[str, Model],
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
        del kwargs
        assert item_data
        v_conf_site_name, data = item_data
        conf_site_name, version = strip_version(v_conf_site_name, separator="+")

        canon_site = canon_sites[v_conf_site_name]

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
            "canon_site": canon_site,
            "version": version,
        }

        defaults = {
            "residues": residues,
        }

        members = extract(key="members")

        ref_ligands = extract(key="reference_ligand_id")

        index_fields = {
            "members": members,
            "reference_ligands": ref_ligands,
        }

        return ProcessedObject(
            model_class=CanonSiteConf,
            fields=fields,
            index_data=index_fields,
            key=conf_site_name,
            versioned_key=v_conf_site_name,
            defaults=defaults,
        )

    @create_objects(depth=1)
    def process_xtalform_site(
        self,
        xtalforms: dict[int | str, MetadataObject],
        canon_sites: dict[str, Model],
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
        del kwargs
        assert item_data
        v_xtalform_site_name, data = item_data
        xtalform_site_name, version = strip_version(v_xtalform_site_name)

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
            "version": version,
        }

        defaults = {
            "lig_chain": lig_chain,
            "residues": residues,
        }

        index_data = {
            "residues": residues,
        }

        return ProcessedObject(
            model_class=XtalformSite,
            fields=fields,
            defaults=defaults,
            key=xtalform_site_name,
            versioned_key=v_xtalform_site_name,
            index_data=index_data,
        )

    @create_objects(depth=6)
    def process_site_observation(
        self,
        experiments: dict[int | str, MetadataObject],
        compounds: dict[int | str, MetadataObject],
        xtalform_sites: dict[str, Model],
        canon_site_confs: dict[int | str, MetadataObject],
        item_data: tuple[str, str, str, int | str, int, str, dict] | None = None,
        # chain: str,
        # ligand: str,
        # version: int,
        # idx: int | str,
        # data: dict,
        validate_files: bool = True,
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
        del kwargs
        assert item_data
        try:
            experiment_id, _, chain, ligand, version, v_idx, data = item_data
        except ValueError:
            # wrong data item
            return None

        extract = functools.partial(
            self._extract,
            data=data,
            section_name="crystals",
            item_name=experiment_id,
            level=logging.WARNING,
        )

        experiment = experiments[experiment_id].instance

        longcode = (
            # f"{experiment.code}_{chain}_{str(ligand)}_{str(version)}_{str(v_idx)}"
            f"{experiment.code}_{chain}_{str(ligand)}_v{str(version)}"
        )
        key = f"{experiment.code}/{chain}/{str(ligand)}"
        v_key = f"{experiment.code}/{chain}/{str(ligand)}/{version}"

        smiles = extract(key="ligand_smiles_string")
        ligand_name = extract(key="ligand_name")

        try:
            compound = compounds[(experiment_id, ligand_name)].instance  # type: ignore[index]
            # I don't understand the error above, I've definitely declared the tuple type
        except KeyError:
            # compound not saved on this round, but if this is not the
            # first upload, experiment and compound may have come from
            # the first.
            try:
                logger.debug('exp: %s, %s', experiment, experiments[experiment_id].new)
                compound = experiment.compounds.get(
                    smiles=experiments[experiment_id].index_data["smiles"]
                )
            except Compound.DoesNotExist:
                # really doensn't exist, can happen
                compound = None
                self.report.log(
                    logging.INFO,
                    f"No compounds for experiment {experiment.code}",
                )
            except MultipleObjectsReturned:
                # based on the way data is presented, I'm fairly
                # certain this cannot happen. but there's nothing in
                # the db to prevent this
                compound = experiment.compounds.filter(smiles=smiles).first()
                self.report.log(
                    logging.WARNING,
                    f"Multiple compounds for experiment {experiment.code}",
                )

        canon_site_conf = canon_site_confs[v_idx].instance
        xtalform_site = xtalform_sites[v_key]

        (  # pylint: disable=unbalanced-tuple-unpacking
            bound_file_t,
            apo_solv_file_t,
            apo_desolv_file_t,
            apo_file_t,
            artefacts_file_t,
            sigmaa_file_t,
            diff_file_t,
            event_file_t,
            ligand_pdb_t,
            ligand_mol_t,
            ligand_smiles_t,
        ) = self.validate_files(
            obj_identifier=experiment_id,
            file_struct=data,
            required=(
                "structure",
                "pdb_apo_solv",
                "pdb_apo_desolv",
                "pdb_apo",
            ),
            recommended=(
                "artefacts",
                "sigmaa_map",  # NB! keys in meta_aligner not yet updated
                "diff_map",  # NB! keys in meta_aligner not yet updated
                "event_map",
                "ligand_pdb",
                "ligand_mol",
                "ligand_smiles",
            ),
            validate_files=validate_files,
        )

        bound_file = bound_file_t[0]
        apo_solv_file = apo_solv_file_t[0]
        apo_desolv_file = apo_desolv_file_t[0]
        apo_file = apo_file_t[0]
        artefacts_file = artefacts_file_t[0]
        sigmaa_file = sigmaa_file_t[0]
        diff_file = diff_file_t[0]
        event_file = event_file_t[0]
        ligand_pdb = ligand_pdb_t[0]
        ligand_mol = ligand_mol_t[0]
        ligand_smiles = ligand_smiles_t[0]

        fields = {
            # Code for this protein (e.g. Mpro_Nterm-x0029_A_501_0)
            "longcode": longcode,
            "version": version,
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
            "sigmaa_file": str(self._get_final_path(sigmaa_file)),
            "diff_file": str(self._get_final_path(diff_file)),
            "event_file": str(self._get_final_path(event_file)),
            "artefacts_file": str(self._get_final_path(artefacts_file)),
            "ligand_pdb": str(self._get_final_path(ligand_pdb)),
            "ligand_mol": str(self._get_final_path(ligand_mol)),
            "ligand_smiles": str(self._get_final_path(ligand_smiles)),
            "pdb_header_file": "currently missing",
        }

        return ProcessedObject(
            model_class=SiteObservation,
            fields=fields,
            defaults=defaults,
            key=key,
            versioned_key=v_key,
        )

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
            msg = "Upload directory missing from uploaded file"
            self.report.log(logging.ERROR, msg)
            # what do you mean unused?!
            raise StopIteration(
                msg
            ) from exc  # pylint: disable=# pylint: disable=protected-access

        with contextlib.suppress(StopIteration):
            upload_dir = next(up_iter)
            self.report.log(logging.WARNING, "Multiple upload directories in archive")
        # now that target name is not included in path, I don't need
        # it here, need it just before creating target object. Also,
        # there's probably no need to throw a fatal here, I can
        # reasonably well deduce it from meta (I think)
        config_it = upload_dir.glob(CONFIG_FILE)
        try:
            config_file = next(config_it)
        except StopIteration as exc:
            msg = f"config file missing from {str(upload_dir)}"
            self.report.log(logging.ERROR, msg)
            raise StopIteration() from exc

        # load necessary files
        config = self._load_yaml(config_file)
        meta = self._load_yaml(Path(upload_dir).joinpath(METADATA_FILE))
        xtalforms_yaml = self._load_yaml(Path(upload_dir).joinpath(XTALFORMS_FILE))

        # this is the last file to load. if any of the files missing, don't continue
        if not any([meta, config, xtalforms_yaml]):
            msg = "Missing files in uploaded data, aborting"
            raise FileNotFoundError(msg)

        # Validate the upload's XCA version information against any MODE-based conditions.
        # An error message is returned if the bundle is not supported.
        if vb_err_msg := _validate_bundle_against_mode(config):
            self.report.log(logging.ERROR, vb_err_msg)
            raise AssertionError(vb_err_msg)

        # Target (very least) is required
        try:
            self.target_name = config["target_name"]
        except KeyError as exc:
            msg = "target_name missing in config file"
            self.report.log(logging.ERROR, msg)
            raise KeyError(msg) from exc

        # project needs to be created before target
        # TODO: original target loader's function get_create_projects
        # seems to handle more cases. adopt or copy
        visit = self.proposal_ref.split()[0]
        self.project, project_created = Project.objects.get_or_create(title=visit)

        self.target, target_created = Target.objects.get_or_create(
            title=self.target_name,
            display_name=self.target_name,
            project=self.project,
        )

        if target_created:
            target_dir = f"{self.target_name}_{self.proposal_ref}"
            # mypy thinks target and target_name are None
            target_dir = sanitize_directory_name(target_dir, self.abs_final_path)  # type: ignore [arg-type]
            self.target.zip_archive = target_dir  # type: ignore [attr-defined]
            self.target.save()  # type: ignore [attr-defined]
        else:
            # NB! using existing field zip_archive to point to the
            # location of the archives, not the archives
            # themselves. The field was unused, and because of the
            # versioned uploads, there's no single archive anymore
            target_dir = str(self.target.zip_archive)  # type: ignore [attr-defined]

        self._final_path = self._final_path.joinpath(target_dir)
        self._abs_final_path = self._abs_final_path.joinpath(target_dir)

        try:
            committer = get_user_model().objects.get(pk=self.user_id)
        except get_user_model().DoesNotExist:
            # add upload as anonymous user
            committer = get_user_model().objects.get(pk=settings.ANONYMOUS_USER)

        # collect top level info
        self.version_number = int(meta["version_number"])
        self.version_dir = meta["version_dir"]
        self.previous_version_dirs = meta["previous_version_dirs"]
        prefix_tooltips = meta.get("code_prefix_tooltips", {})

        # TODO: is it here where I can figure out if this has already been uploaded?
        if self._is_already_uploaded(target_created, project_created):
            # remove uploaded file
            Path(self.bundle_path).unlink()
            msg = f"{self.bundle_name} already uploaded"
            self.report.log(logging.ERROR, msg)
            raise FileExistsError(msg)

        if project_created and self.project.title in settings.PUBLIC_TAS_LIST:  # type: ignore[attr-defined]
            assert self.project
            self.project.open_to_public = True
            self.project.save()

        assert self.target

        # check transformation matrix files
        (  # pylint: disable=unbalanced-tuple-unpacking
            trans_neighbourhood,
            trans_conf_site,
            trans_ref_struct,
        ) = self.validate_files(
            obj_identifier="trans_matrices",
            # since the paths are given if file as strings, I think I
            # can get away with compiling them as strings here
            file_struct={
                TRANS_NEIGHBOURHOOD: f"{self.version_dir}/{TRANS_NEIGHBOURHOOD}",
                TRANS_CONF_SITE: f"{self.version_dir}/{TRANS_CONF_SITE}",
                TRANS_REF_STRUCT: f"{self.version_dir}/{TRANS_REF_STRUCT}",
            },
            required=(TRANS_NEIGHBOURHOOD, TRANS_CONF_SITE, TRANS_REF_STRUCT),
        )

        trans_neighbourhood = trans_neighbourhood[0]
        trans_conf_site = trans_conf_site[0]
        trans_ref_struct = trans_ref_struct[0]

        self.experiment_upload.project = self.project
        self.experiment_upload.target = self.target
        self.experiment_upload.committer = committer
        self.experiment_upload.neighbourhood_transforms = str(
            self._get_final_path(trans_neighbourhood)
        )
        self.experiment_upload.conformer_site_transforms = str(
            self._get_final_path(trans_conf_site)
        )
        self.experiment_upload.reference_structure_transforms = str(
            self._get_final_path(trans_ref_struct)
        )
        self.experiment_upload.upload_data_dir = self.version_dir
        self.experiment_upload.upload_version = self.version_number
        self.experiment_upload.save()

        (  # pylint: disable=unbalanced-tuple-unpacking
            assemblies,
            xtalform_assemblies,
        ) = self._get_yaml_blocks(
            yaml_data=xtalforms_yaml,
            # blocks=("assemblies", "xtalforms"),
            blocks=("assemblies", "crystalforms"),
        )

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
                "crystalforms",
                "canon_sites",
                "conformer_sites",
                "xtalform_sites",
            ),
        )

        experiment_objects = self.process_experiment(
            yaml_data=crystals, prefix_tooltips=prefix_tooltips
        )
        compound_objects = self.process_compound(
            yaml_data=crystals, experiments=experiment_objects
        )

        # save components manytomany to experiment
        # TODO: is it 1:1 relationship? looking at the meta_align it
        # seems to be, but why the m2m then?
        for (
            comp_code,  # it's tuple here, (epx_name, ligand_key)
            comp_meta,
        ) in compound_objects.items():  # pylint: disable=no-member
            experiment = experiment_objects[comp_code[0]].instance
            experiment.compounds.add(comp_meta.instance)
            comp_meta.instance.project_id.add(self.experiment_upload.project)

        xtalform_objects = self.process_xtalform(yaml_data=xtalforms)
        self._enumerate_objects(xtalform_objects, "xtalform_num")

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
        self._enumerate_objects(quat_assembly_objects, "assembly_num")

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
        # TODO: ref_conf_site is with version, object's key isn't
        for val in canon_site_objects.values():  # pylint: disable=no-member
            val.instance.ref_conf_site = canon_site_conf_objects[
                val.index_data["reference_conformer_site_id"]
            ].instance
            val.instance.save()

        # canon site instances are now complete
        # still missing fk to site_observation in canon_site_conf

        # reindex xtalform site to grab for site observation
        xtalform_site_by_tag = {}
        for val in xtalform_sites_objects.values():  # pylint: disable=no-member
            for k in val.index_data["residues"]:
                xtalform_site_by_tag[k] = val.instance

        site_observation_objects = self.process_site_observation(
            yaml_data=crystals,
            experiments=experiment_objects,
            compounds=compound_objects,
            xtalform_sites=xtalform_site_by_tag,
            canon_site_confs=canon_site_conf_objects,
        )

        values = ["experiment"]
        # fmt: off
        qs = SiteObservation.objects.filter(
                experiment__experiment_upload__target=self.target,
                code__isnull=True,
            ).values(
                *values,
            ).order_by(
                *values,
            ).annotate(
                obvs=ArrayAgg("id"),
            ).values_list("obvs", flat=True)
        # fmt: on

        for elem in qs:
            # fmt: off
            subgroups = SiteObservation.objects.filter(
                pk__in=elem,
            ).order_by(
                "canon_site_conf__canon_site",
            ).annotate(
                sites=Count("canon_site_conf__canon_site"),
                obvs=ArrayAgg('id'),
            ).order_by(
                "-sites",
            ).values_list("obvs", flat=True)
            # fmt: on

            suffix = alphanumerator()
            for sub in subgroups:
                # objects in this group should be named with same scheme
                so_group = SiteObservation.objects.filter(pk__in=sub)

                # memo to self: there used to be some code to test the
                # position of the iterator in existing entries. This
                # was because it was assumed, that when adding v2
                # uploads, it can bring along new observations under
                # existing experiment. Following discussions with
                # Conor, it seems that this will not be the case. But
                # should it agin be, this code was deleted on
                # 2024-03-04, if you need to check

                for so in so_group.filter(code__isnull=True):
                    logger.debug("processing so: %s", so.longcode)
                    if so.experiment.type == 1:
                        # manual. code is pdb code
                        code = f"{so.experiment.code}-{next(suffix)}"
                        # NB! at the time of writing this piece of
                        # code, I haven't seen an example of the data
                        # so I only have a very vague idea how this is
                        # going to work. The way I understand it now,
                        # they cannot belong to separate groups so
                        # there's no need for different iterators. But
                        # could be I need to split them up
                    else:
                        # model building. generate code
                        code_prefix = experiment_objects[so.experiment.code].index_data[
                            "code_prefix"
                        ]
                        # iter_pos = next(suffix)
                        # code = f"{code_prefix}{so.experiment.code.split('-')[1]}{iter_pos}"
                        # code = f"{code_prefix}{so.experiment.code.split('-')[1]}{next(suffix)}"
                        code = f"{code_prefix}{so.experiment.code.split('-x')[1]}{next(suffix)}"

                        # test uniqueness for target
                        # TODO: this should ideally be solved by db engine, before
                        # rushing to write the trigger, have think about the
                        # loader concurrency situations
                        code_qs = SiteObservation.objects.filter(
                            experiment__experiment_upload__target=self.target,
                            code=code,
                        )
                        # if code exists and the experiment is new
                        logger.debug(
                            'checking code uniq: %s, %s', code, so.experiment.status
                        )
                        if code_qs.exists() and so.experiment.status == 0:
                            msg = (
                                f"short code {code} already exists for this target; "
                                + "specify a code_prefix to resolve this conflict"
                            )
                            self.report.log(logging.ERROR, msg)

                    so.code = code
                    so.save()

        # site_observations_versioned = {}
        # for val in site_observation_objects.values():  # pylint: disable=no-member
        #     site_observations_versioned[val.versioned_key] = val.instance

        # final remaining fk, attach reference site observation to canon_site_conf
        for val in canon_site_conf_objects.values():  # pylint: disable=no-member
            val.instance.ref_site_observation = site_observation_objects[
                val.index_data["reference_ligands"]
            ].instance
            logger.debug("attaching canon_site_conf: %r", val.instance)
            logger.debug(
                "attaching canon_site_conf: %r",
                site_observation_objects[
                    val.index_data["reference_ligands"]
                ].instance.longcode,
            )
            val.instance.save()

        logger.debug("data read and processed, adding tags")

        # tag site observations
        cat_canon = TagCategory.objects.get(category="CanonSites")
        # sort canon sites by number of observations
        # fmt: off
        canon_sort_qs = CanonSite.objects.filter(
            pk__in=[k.instance.pk for k in canon_site_objects.values() ], # pylint: disable=no-member
        ).annotate(
            # obvs=Count("canonsiteconf_set__siteobservation_set", default=0),
            obvs=Count("canonsiteconf__siteobservation", default=0),
        ).order_by("-obvs", "name")
        # ordering by name is not strictly necessary, but
        # makes the sorting consistent

        # fmt: on

        logger.debug('canon_site_order')
        for site in canon_sort_qs:
            logger.debug('%s: %s', site.name, site.obvs)

        _canon_site_objects = {}
        for site in canon_sort_qs:
            key = f"{site.name}+{site.version}"
            _canon_site_objects[key] = canon_site_objects[
                key
            ]  # pylint: disable=no-member

        self._enumerate_objects(_canon_site_objects, "canon_site_num")
        for val in _canon_site_objects.values():  # pylint: disable=no-member
            prefix = val.instance.canon_site_num
            # tag = canon_name_tag_map.get(val.versioned_key, "UNDEFINED")
            so_list = SiteObservation.objects.filter(
                canon_site_conf__canon_site=val.instance
            )
            tag = val.versioned_key
            try:
                short_tag = val.versioned_key.split('-')[1][1:]
                main_obvs = val.instance.ref_conf_site.ref_site_observation
                code_prefix = experiment_objects[main_obvs.experiment.code].index_data[
                    "code_prefix"
                ]
                short_tag = f"{code_prefix}{short_tag}"
            except IndexError:
                short_tag = tag

            self._tag_observations(
                tag,
                prefix,
                category=cat_canon,
                site_observations=so_list,
                short_tag=short_tag,
            )

        logger.debug("canon_site objects tagged")

        numerators = {}
        cat_conf = TagCategory.objects.get(category="ConformerSites")
        for val in canon_site_conf_objects.values():  # pylint:
            # disable=no-member problem introduced with the sorting of
            # canon sites (issue 1498). objects somehow go out of sync
            val.instance.refresh_from_db()
            if val.instance.canon_site.canon_site_num not in numerators.keys():
                numerators[val.instance.canon_site.canon_site_num] = alphanumerator()
            prefix = (
                f"{val.instance.canon_site.canon_site_num}"
                + f"{next(numerators[val.instance.canon_site.canon_site_num])}"
            )
            so_list = [
                site_observation_objects[k].instance for k in val.index_data["members"]
            ]
            # tag = val.instance.name.split('+')[0]
            tag = val.instance.name
            try:
                short_tag = val.instance.name.split('-')[1][1:]
                main_obvs = val.instance.ref_site_observation
                code_prefix = experiment_objects[main_obvs.experiment.code].index_data[
                    "code_prefix"
                ]
                short_tag = f"{code_prefix}{short_tag}"
            except IndexError:
                short_tag = tag

            self._tag_observations(
                tag,
                prefix,
                category=cat_conf,
                site_observations=so_list,
                hidden=True,
                short_tag=short_tag,
            )

        logger.debug("conf_site objects tagged")

        cat_quat = TagCategory.objects.get(category="Quatassemblies")
        for val in quat_assembly_objects.values():  # pylint: disable=no-member
            prefix = f"A{val.instance.assembly_num}"
            tag = val.instance.name
            so_list = SiteObservation.objects.filter(
                xtalform_site__xtalform__in=XtalformQuatAssembly.objects.filter(
                    quat_assembly=val.instance
                ).values("xtalform")
            )
            self._tag_observations(
                tag, prefix, category=cat_quat, site_observations=so_list
            )

        logger.debug("quat_assembly objects tagged")

        cat_xtal = TagCategory.objects.get(category="Crystalforms")
        for val in xtalform_objects.values():  # pylint: disable=no-member
            prefix = f"F{val.instance.xtalform_num}"
            so_list = SiteObservation.objects.filter(
                xtalform_site__xtalform=val.instance
            )
            tag = val.instance.name

            self._tag_observations(
                tag, prefix, category=cat_xtal, site_observations=so_list
            )

        logger.debug("xtalform objects tagged")

        # enumerate xtalform_sites. a bit trickier than others because
        # requires alphabetic enumeration starting from the letter of
        # the chain and following from there

        # sort the dictionary
        # fmt: off
        xtls_sort_qs = XtalformSite.objects.filter(
            pk__in=[k.instance.pk for k in xtalform_sites_objects.values() ], # pylint: disable=no-member
        ).annotate(
            obvs=Count("canon_site__canonsiteconf__siteobservation", default=0),
        ).order_by("-obvs", "xtalform_site_id")
        # ordering by xtalform_site_id is not strictly necessary, but
        # makes the sorting consistent

        # fmt: on

        _xtalform_sites_objects = {}
        for xtl in xtls_sort_qs:
            key = f"{xtl.xtalform_site_id}/{xtl.version}"
            _xtalform_sites_objects[key] = xtalform_sites_objects[
                key
            ]  # pylint: disable=no-member

        if self.version_number == 1:
            # first upload, use the chain letter
            xtnum = alphanumerator(
                start_from=xtls_sort_qs[0].lig_chain.lower(), drop_first=False
            )
        else:
            # subsequent upload, just use the latest letter as starting point
            # fmt: off
            last_xtsite = XtalformSite.objects.filter(
                    pk__in=[
                        k.instance.pk
                        for k in _xtalform_sites_objects.values()  # pylint: disable=no-member
                    ]
                ).order_by(
                    "-xtalform_site_num"
                )[0].xtalform_site_num
            # fmt: on
            xtnum = alphanumerator(start_from=last_xtsite)

            # this should be rare, as Frank said, all crystal-related
            # issues should be resolved by the time of the first
            # upload. In fact, I'll mark this momentous occasion here:
            logger.warning("New XtalformSite objects added in subsequent uploads")

        for val in _xtalform_sites_objects.values():  # pylint: disable=no-member
            if not val.instance.xtalform_site_num:
                val.instance.xtalform_site_num = next(xtnum)
                val.instance.save()

        cat_xtalsite = TagCategory.objects.get(category="CrystalformSites")
        for val in _xtalform_sites_objects.values():  # pylint: disable=no-member
            prefix = (
                f"F{val.instance.xtalform.xtalform_num}"
                + f"{val.instance.xtalform_site_num}"
            )
            so_list = [
                site_observation_objects[k].instance for k in val.index_data["residues"]
            ]
            tag = val.versioned_key
            try:
                # remove protein name and 'x'
                short_tag = val.instance.xtalform_site_id.split('-')[1][1:]
                main_obvs = val.instance.canon_site.ref_conf_site.ref_site_observation
                code_prefix = experiment_objects[main_obvs.experiment.code].index_data[
                    "code_prefix"
                ]
                short_tag = f"{code_prefix}{short_tag}"
            except IndexError:
                short_tag = tag

            self._tag_observations(
                tag,
                prefix,
                category=cat_xtalsite,
                site_observations=so_list,
                hidden=True,
                short_tag=short_tag,
            )

        logger.debug("xtalform_sites objects tagged")

        self._generate_poses()

        # tag all new observations, so that the curator can find and
        # re-pose them
        self._tag_observations(
            self.version_dir,
            "",
            TagCategory.objects.get(category="Other"),
            [
                k.instance
                for k in site_observation_objects.values()  # pylint: disable=no-member
                if k.new
            ],
            clean_ids=False,
        )

    def _load_yaml(self, yaml_file: Path) -> dict:
        contents = {}
        try:
            with open(yaml_file, "r", encoding="utf-8") as file:
                contents = yaml.safe_load(file)
        except FileNotFoundError:
            self.report.log(
                logging.ERROR, f"File {yaml_file.name} not found in data archive"
            )

        return contents

    # TODOL error handling. what's the correct response when
    # something's missing? push through and compile report?
    def _get_yaml_blocks(self, yaml_data: dict, blocks: Iterable) -> list[dict]:
        error_text = "'{}' section missing in input file"
        result = []
        for block in blocks:
            try:
                result.append(yaml_data[block])
            except KeyError:
                msg = error_text.format(block)
                self.report.log(logging.ERROR, msg)

        return result

    def _extract(
        self,
        data: dict,
        key: str | int,
        section_name: str,
        item_name: str,
        level: int = logging.ERROR,
        return_type: type = str,
    ) -> Any:
        try:
            result = data[key]
        except KeyError as exc:
            result = "" if level == logging.INFO else "missing"
            if return_type == list:
                result = [result]

            self.report.log(
                level, f"{exc} missing from {section_name}: {item_name} section"
            )

        return result

    def _generate_poses(self):
        values = ["canon_site_conf__canon_site", "cmpd"]
        # fmt: off
        pose_groups = SiteObservation.filter_manager.by_target(
            self.target,
        ).exclude(
            canon_site_conf__canon_site__isnull=True,
        ).exclude(
            cmpd__isnull=True,
        ).values(
            *values
        ).order_by(
            "canon_site_conf__canon_site",
        ).annotate(
            obvs=ArrayAgg('id'),
        ).values_list("obvs", flat=True)
        # fmt: on

        for group in pose_groups:
            pose_items = SiteObservation.objects.filter(pk__in=group)
            sample = pose_items.first()
            # check for existing group
            try:
                pose = Pose.objects.get(
                    canon_site=sample.canon_site_conf.canon_site,
                    compound=sample.cmpd,
                )
            except Pose.DoesNotExist:
                # create new, add random observation as main
                pose = Pose(
                    canon_site=sample.canon_site_conf.canon_site,
                    compound=sample.cmpd,
                    main_site_observation=sample,
                    display_name=sample.code,
                )
                pose.save()
            except MultipleObjectsReturned:
                # must be a follow-up upload. create new pose, but
                # only add observatons that are not yet assigned (if
                # these exist)
                pose_items = pose_items.filter(pose__isnull=True)
                if pose_items.exists():
                    sample = pose_items.first()
                    pose = Pose(
                        canon_site=sample.canon_site_conf.canon_site,
                        compound=sample.cmpd,
                        main_site_observation=sample,
                        display_name=sample.code,
                    )
                    pose.save()
                else:
                    # I don't know if this can happen but this (due to
                    # other bugs) is what allowed me to find this
                    # error. Make a note in the logs.
                    logger.warning("No observations left to assign to pose")

            # finally add observations to the (new or existing) pose
            for obvs in pose_items:
                obvs.pose = pose
                obvs.save()

    def _tag_observations(
        self,
        tag: str,
        prefix: str,
        category: TagCategory,
        site_observations: list,
        hidden: bool = False,
        short_tag: str | None = None,
        clean_ids: bool = True,
    ) -> None:
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
            assert self.target
            so_group = SiteObservationGroup(target=self.target)
            so_group.save()
        except MultipleObjectsReturned:
            SiteObservationGroup.objects.filter(
                target=self.target, description=tag
            ).delete()
            assert self.target
            so_group = SiteObservationGroup(target=self.target)
            so_group.save()

        name = f"{prefix} - {tag}" if prefix else tag
        tag = tag if short_tag is None else short_tag
        short_name = name if short_tag is None else f"{prefix} - {short_tag}"

        if clean_ids:
            tag = clean_object_id(tag)
            name = clean_object_id(name)
            short_name = clean_object_id(short_name)

        try:
            so_tag = SiteObservationTag.objects.get(
                upload_name=name, target=self.target
            )
            # Tag already exists
            # Apart from the new mol_group and molecules, we shouldn't be
            # changing anything.
            so_tag.mol_group = so_group
        except SiteObservationTag.DoesNotExist:
            so_tag = SiteObservationTag(
                tag=tag,
                tag_prefix=prefix,
                upload_name=name,
                category=category,
                target=self.target,
                mol_group=so_group,
                hidden=hidden,
                short_tag=short_name,
            )

        so_tag.save()

        so_group.site_observation.add(*site_observations)
        so_tag.site_observations.add(*site_observations)

    def _is_already_uploaded(self, target_created, project_created):
        if target_created or project_created:
            return False
        else:
            uploaded = ExperimentUpload.objects.filter(
                target=self.target,
                project=self.project,
            ).values_list("upload_data_dir", flat=True)

            return self.version_dir in uploaded

    def _get_final_path(self, path: str | None) -> Path | None:
        """Update relative path to final storage path

        NB! this returns a relative path that can be used in queries
        not absoulte one. This is used to populate location fields in
        database tables.
        """
        try:
            return self.final_path.joinpath(path)  # type: ignore[arg-type]
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
    del contact_email
    with TemporaryDirectory(dir=settings.MEDIA_ROOT) as tempdir:
        target_loader = TargetLoader(
            data_bundle, proposal_ref, tempdir, user_id=user_id, task=task
        )

        # Decompression can take some time, so we want to report progress
        bundle_filename = os.path.basename(data_bundle)
        target_loader.report.log(logging.INFO, f"Decompressing '{bundle_filename}'")

        try:
            # archive is first extracted to temporary dir and moved later
            with tarfile.open(target_loader.bundle_path, "r") as archive:
                msg = f"Extracting bundle: {data_bundle}"
                logger.info("%s%s", target_loader.report.task_id, msg)
                archive.extractall(target_loader.raw_data)
                msg = f"Data extraction complete: {data_bundle}"
                logger.info("%s%s", target_loader.report.task_id, msg)
        except Exception as exc:
            # Handle _any_ underlying problem with the file.
            logger.error('Got an exception opening the file: %s', str(exc))
            target_loader.report.log(
                logging.ERROR,
                f"Decompression of '{bundle_filename}' has failed. Is it a Target Experiment file?",
            )
            target_loader.report.final(
                f"Failed to decompress '{target_loader.data_bundle}'", success=False
            )
            return

        target_loader.report.log(logging.INFO, f"Decompressed '{bundle_filename}'")

        try:
            with transaction.atomic():
                target_loader.process_bundle()
                if target_loader.report.failed:
                    # need to trigger transaction failure
                    raise IntegrityError(
                        f"Uploading {target_loader.data_bundle} failed"
                    )
        except Exception as exc:
            # Handle _any_ underlying problem.
            # These are errors processing the data, which we handle gracefully.
            # The task should _always_ end successfully.
            # Any problem with the underlying data is transmitted in the report.
            logger.error(exc, exc_info=True)
            target_loader.report.final(
                f"Failed to process '{target_loader.data_bundle}'", success=False
            )
            return
        else:
            _move_and_save_target_experiment(target_loader)


def _move_and_save_target_experiment(target_loader):
    # Move the uploaded file to its final location
    try:
        target_loader.abs_final_path.mkdir(parents=True)
    except FileExistsError:
        # subsequent upload, directory already exists
        pass

    shutil.move(
        str(target_loader.raw_data.joinpath(target_loader.version_dir)),
        str(target_loader.abs_final_path),
    )
    Path(target_loader.bundle_path).rename(
        target_loader.abs_final_path.joinpath(target_loader.data_bundle)
    )

    set_directory_permissions(target_loader.abs_final_path, 0o755)

    target_loader.report.final(f"{target_loader.data_bundle} uploaded successfully")
    target_loader.experiment_upload.message = target_loader.report.json()
    target_loader.experiment_upload.save()


# when moving to class:
# def import_compound_identifiers(self):
def import_compound_identifiers(filename: str):
    try:
        # df = pd.read_csv(CUSTOM_IDENTIFIER_FILE)
        df = pd.read_csv(filename)
    except UnicodeDecodeError:
        logger.debug("file import error")
        # self.report.log(
        #     logging.ERROR,
        #     f"Error reading {CUSTOM_IDENTIFIER_FILE}, unexpected format",
        # )
        return

    identifiers_from_file = set(
        [k for k in df.columns if k not in CUSTOM_IDENTIFIER_COLUMNS]
    )

    # I think this is a bad idea, but it was explicitly in the spec
    identifier_types = set(
        CompoundIdentifierType.objects.values_list("name", flat=True)
    )
    new_identifiers = identifiers_from_file.difference(identifier_types)
    for identifier in new_identifiers:
        CompoundIdentifierType(name=identifier).save()

    # you'd think I could supply the compounds processed, but I need a queryset..
    compounds = Compound.objects.annotate(
        exp_code=F("experiment__code"),
    ).filter(
        experiment__code__in=df["xtal"],
    )

    identifiers = CompoundIdentifierType.objects.all()

    for idf in identifiers_from_file:
        identifier = identifiers.get(name=idf)

        for _, row in df.loc[
            df[idf].notna(), list(CUSTOM_IDENTIFIER_COLUMNS) + [idf]
        ].iterrows():
            exp_code, ligand_name, compound_code, name = row
            compound = compounds.get(exp_code=exp_code, ligand_name=ligand_name)

            if compound.compound_code == compound_code:
                CompoundIdentifier(
                    type=identifier,
                    compound=compound,
                    name=name,
                ).save()
            else:
                pass
                # self.report.log(
                #     logging.ERROR,
                #     f"compound_code changed for {exp_code}, {ligand_name}",
                # )

            # memo to self: I tried using bulk_create here but it kept
            # failing with a rather cryptic error message. worth
            # revisiting when django has been upgraded
