import contextlib
import copy
import functools
import hashlib
import logging
import math
import os
import re
import shutil
import sqlite3
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import timedelta
from enum import Enum
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Iterable, List, Optional, Tuple, TypeVar

import pandas as pd
import yaml
from celery import Task
from dateutil.parser import parse
from dateutil.parser._parser import ParserError  # type: ignore [import-untyped]
from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.postgres.aggregates import ArrayAgg
from django.core.exceptions import MultipleObjectsReturned
from django.db import IntegrityError, transaction
from django.db.models import Count, F, Model
from django.db.models.base import ModelBase
from django.utils import timezone
from rdkit import Chem

from api.utils import deployment_mode_is_production
from fragalysis.settings import TARGET_LOADER_MEDIA_DIRECTORY
from viewer.models import (
    AtomCoordinates,
    CanonSite,
    CanonSiteConf,
    Compound,
    CompoundIdentifier,
    CompoundIdentifierType,
    Experiment,
    ExperimentStatusType,
    ExperimentUpload,
    Pose,
    Project,
    QualityStatusType,
    QuatAssembly,
    SiteObservation,
    SiteObservationComputedSiteObservation,
    SiteObservationQualityStatus,
    TagCategory,
    Target,
    Xtalform,
    XtalformQuatAssembly,
    XtalformSite,
)
from viewer.utils import alphanumerator, longcode_from_tag, sanitize_directory_name

from .tags import TagManager

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


def _get_xca_tag(yaml_content: dict[str, Any]) -> tuple[list[str], str]:
    # Initial concern - the loader's git information.
    # It must not be 'dirty' and must have a valid 'tag'.
    xca_git_info_key = "xca_git_info"
    base_error_msg = "Stack is in PRODUCTION mode - and"
    try:
        xca_git_info = yaml_content[xca_git_info_key]
    except KeyError as exc:
        raise ValueError(
            f"{base_error_msg} '{xca_git_info_key}' is a required configuration property"
        ) from exc

    logger.info("%s: %s", xca_git_info_key, xca_git_info)

    if "dirty" not in xca_git_info:
        raise ValueError(
            f"{base_error_msg} '{xca_git_info_key}' has no 'dirty' property"
        )
    if xca_git_info["dirty"]:
        raise ValueError(f"{base_error_msg} '{xca_git_info_key}->dirty' must be False")

    if "tag" not in xca_git_info:
        raise ValueError(f"{base_error_msg} '{xca_git_info_key}' has no 'tag' property")
    xca_version_tag: str = str(xca_git_info["tag"])
    tag_parts: List[str] = xca_version_tag.split(".")

    return tag_parts, xca_version_tag


def _check_xca_tag(
    yaml_content: dict[str, Any],
    min_xca_tag: str,
) -> bool:
    if not deployment_mode_is_production():
        # We're not in production mode - no bundle checks
        return True

    logger.debug("Checking XCA tag")
    try:
        tag_parts, _ = _get_xca_tag(yaml_content)
    except ValueError as exc:
        raise ValueError(exc.args[0]) from exc

    min_tag = [int(k) for k in min_xca_tag.split(".")]

    logger.debug("Given tag: %s", tag_parts)
    logger.debug("Settings min tag: %s", min_xca_tag)

    try:
        xca_tag: list[int] = [int(k) for k in tag_parts[:2]]
    except ValueError as exc:
        raise ValueError(f"{'.'.join(tag_parts)} is not a valid tag") from exc

    for k, v in zip(xca_tag, min_tag):
        logger.debug("XCA     tag checker: k=%s, v=%s", k, v)
        if k > v:
            return False

    return True


def _validate_bundle_against_mode(config_yaml: Dict[str, Any]) -> Optional[str]:
    """Inspects the meta to ensure it is supported by the MODE this stack is in.
    Mode is (typically) one of DEVELOPER or PRODUCTION.
    """
    assert config_yaml
    if not deployment_mode_is_production():
        # We're not in production mode - no bundle checks
        return None

    # PRODUCTION mode (strict)
    # duplicated to produce the error message
    # TODO: better refactor
    xca_git_info_key = "xca_git_info"
    base_error_msg = "Stack is in PRODUCTION mode - and"
    try:
        tag_parts, xca_version_tag = _get_xca_tag(config_yaml)
    except ValueError as exc:
        return exc.args[0]

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


def strip_exp_code(code: str) -> str:
    try:
        return re.split(r"-\w{1}", code)[1]
    except IndexError as exc:
        raise ValueError(f"Non-standard experiment code {code}") from exc


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
        logger.debug("+wrapper_create_objects")
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
                # NB!! this is a hack to prevent overwriting
                # experiment_upload value in experiment objects. this
                # only works because no other object passing through
                # here has a field called 'experiment_upload'
                # Alternative: allow NULL and fill later. would that be better?
                instance_data.defaults.pop("experiment_upload", None)

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


def split_version(version_number: str) -> tuple[int, int]:
    splits = version_number.split('.')

    if len(splits) != 2:
        raise ValueError("Unrecognised data format, should be <major>.<minor>")

    try:
        major = int(splits[0])
        minor = int(splits[1])
    except ValueError as exc:
        raise ValueError(f"Non-numeric version number: {version_number}") from exc

    return major, minor


def validate_data_version(
    major: int,
    minor: int,
    o_major: int | None = None,
    o_minor: int | None = None,
    target_name: str | None = None,
    project_name: str | None = None,
) -> Tuple[bool, str]:
    logger.debug('major: %s; minor: %s', major, minor)
    logger.debug('o_major: %s; o_minor: %s', o_major, o_minor)

    s_major, s_minor = [int(k) for k in settings.XCA_DATA_FORMAT_VERSION.split('.')]
    logger.debug('s_major: %s; s_minor: %s', s_major, s_minor)

    if major != s_major:
        return (
            False,
            f"Data major version mismatch: '{s_major}' "
            + f"expected, '{major}' uploaded",
        )

    # alternatively, if target- and project name are given (likely pre-upload check):
    if target_name and project_name and not o_major and not o_minor:
        previous_uploads = ExperimentUpload.objects.filter(
            target__title=target_name,
            project__title=project_name,
        )
        if previous_uploads.exists():
            last_upload = previous_uploads.order_by("upload_version").last()
            o_major = last_upload.data_version_major
            o_minor = last_upload.data_version_minor

    if o_major and o_major < major:
        return False, (
            f"Incoming data major version '{major}' does not match previous upload: "
            + f"'{o_major}'. Please delete the target and prepare new upload"
        )

    if minor != s_minor:
        return (
            True,
            f"Data minor version mismatch: {settings.XCA_DATA_FORMAT_VERSION} "
            + f"expected, {major}.{minor} uploaded",
        )

    if o_minor and o_minor < minor:
        return True, (
            f"Incoming data minor version '{minor}' does not match previous upload: "
            + f"'{o_minor}'"
        )

    # absolutely nothing went wrong
    return True, ""


def validate_upload_version(
    upload_version: int,
    previous_version: int | None = None,
    target_name: str | None = None,
    project_name: str | None = None,
) -> Tuple[bool, str]:
    if not previous_version:
        previous_uploads = ExperimentUpload.objects.filter(
            target__title=target_name,
            project__title=project_name,
        )
        if previous_uploads.exists():
            previous_version = (
                previous_uploads.order_by("upload_version").last().upload_version
            )
        else:
            previous_version = 0

    if previous_version + 1 != upload_version:  # type: ignore [operator]
        return False, (
            f"Upload version {upload_version} is not the expected next version."  # type: ignore [operator]
            f" The next version should be {previous_version + 1}."  # type: ignore [operator]
        )

    # absolutely nothing went wrong
    return True, ""


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
        self.excluded_crystals: list[str] = []

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
            logger.debug("missing file: %s", file_path)
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

        try:
            status = ExperimentStatusType.objects.get(status=dstatus)
        except ExperimentStatusType.DoesNotExist:
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
            # overwrites exp upload in old instances, there's a hack
            # in create_objects method to prevent that
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
            logger.debug("wrong data item")
            return None

        # more validation
        if (
            item_data[1] != "crystallographic_files"
            or item_data[2] != "ligand_cif"
            or item_data[3] != "ligands"
        ):
            logger.debug(
                "wrong data item: %s; %s; %s",
                item_data[1],
                item_data[2],
                item_data[3],
            )
            return None

        smiles = data.get("smiles", None)
        compound_code = data.get("compound_code", None)

        if smiles is None and compound_code is None:
            # gotta have at least something
            logger.debug("no smiles and compound_code")
            return None

        modeled_smiles_soakdb = data.get("modeled_smiles_soakdb", None)
        modeled_smiles_canon = data.get("modeled_smiles_canon", None)
        soaked_smiles_soakdb = data.get("soaked_smiles_soakdb", None)
        soaked_smiles_canon = data.get("soaked_smiles_canon", None)

        inchi_key = ""
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol:
            Chem.RemoveStereochemistry(mol)
            inchi_key = Chem.inchi.MolToInchiKey(mol)

        defaults = {
            "smiles": smiles,
            "inchi_key": inchi_key,
            "compound_code": compound_code,
            "ligand_name": ligand_key,
            "modeled_smiles_soakdb": modeled_smiles_soakdb,
            "modeled_smiles_canon": modeled_smiles_canon,
            "soaked_smiles_soakdb": soaked_smiles_soakdb,
            "soaked_smiles_canon": soaked_smiles_canon,
        }

        fields = {}

        if not experiments[experiment_name].new:
            logger.debug("old experiment: %s", experiment_name)

            try:
                exp = Experiment.filter_manager.by_target(
                    self.target,
                ).get(
                    code=experiment_name,
                    status__isnull=False,
                )
                logger.debug("found old experiment: %s", experiment_name)
                exp_compounds = exp.compounds.all()

                try:
                    cmpd = exp_compounds.get(smiles=smiles)
                    fields = {"id": cmpd.id}
                    logger.debug("found compounds for old experiment: %s", cmpd.pk)
                except Compound.DoesNotExist:
                    logger.debug("did not find compounds for old experiment")
            except Experiment.DoesNotExist:
                logger.debug(
                    "did not find old experiment: %s, it's likely from soqkdb",
                    experiment_name,
                )

        return ProcessedObject(
            model_class=Compound,
            fields=fields,
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
                    ligand_name=ligand_name,
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
            ligand_sdf_t,
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
                "ligand_sdf",
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
        ligand_sdf = ligand_sdf_t[0]

        fields = {
            # Code for this protein (e.g. Mpro_Nterm-x0029_A_501_0)
            # "longcode": longcode,
            "version": version,
            "experiment": experiment,
            "cmpd": compound,
            "xtalform_site": xtalform_site,
            "canon_site_conf": canon_site_conf,
            # "smiles": smiles,
            "seq_id": ligand,
            "chain_id": chain,
        }

        # smiles removed from check fields aand removed to defaults as
        # part of 1670
        # longcode removed as part of 1672, because broke superseding

        defaults = {
            "longcode": longcode,
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
            "ligand_sdf": str(self._get_final_path(ligand_sdf)),
            "pdb_header_file": None,
            "smiles": smiles,
        }

        mol = None
        if ligand_mol:
            molpath = Path(settings.MEDIA_ROOT).joinpath(
                self.raw_data,
                Path(ligand_mol),
            )
            if molpath.exists():
                mol = Chem.MolFromMolFile(str(molpath))

        return ProcessedObject(
            model_class=SiteObservation,
            fields=fields,
            defaults=defaults,
            key=key,
            versioned_key=v_key,
            index_data={'mol': mol},
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
        if vb_err_msg := _validate_bundle_against_mode(meta):
            self.report.log(logging.ERROR, vb_err_msg)
            raise AssertionError(vb_err_msg)

        # check that soakdb output exists
        soakdb_it = Path(upload_dir).joinpath("extra_files").glob("*.sqlite")
        try:
            soakdb_path = next(soakdb_it)
        except StopIteration as exc:
            soakdb_path = None
            msg = f"SoakDB file missing from {Path(upload_dir).joinpath('extra_files')}"
            # does not exist but maybe that's OK?
            try:
                tag_ok = _check_xca_tag(meta, settings.XCA_MIN_VERSION)
                # tag_ok means that the data was compiled with XCA
                # version greater than that specfiied in settings
            except ValueError as e:
                # unsuitable tag
                tag_ok = False
                self.report.log(logging.ERROR, e.args[0])
            if tag_ok:
                self.report.log(logging.WARNING, msg)
            else:
                # version check failed
                self.report.log(logging.ERROR, msg)
                raise StopIteration() from exc

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
            project=self.project,
        )

        if target_created:
            target_dir = f"{self.target_name}_{self.proposal_ref}"
            # mypy thinks target and target_name are None
            target_dir = sanitize_directory_name(target_dir, self.abs_final_path)  # type: ignore [arg-type]
            self.target.zip_archive = target_dir  # type: ignore [attr-defined]
            self.target.display_name = self.target_name  # type: ignore [attr-defined]
            self.target.short_name = self.target_name  # type: ignore [attr-defined]
            self.target.save()  # type: ignore [attr-defined]
        else:
            # NB! using existing field zip_archive to point to the
            # location of the archives, not the archives
            # themselves. The field was unused, and because of the
            # versioned uploads, there's no single archive anymore
            target_dir = str(self.target.zip_archive)  # type: ignore [attr-defined]
            # don't need this on first upload
            self.excluded_crystals = meta.get("excluded_crystals", [])

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
        data_format_version = str(meta["data_format_version"])

        try:
            major, minor = split_version(data_format_version)
        except ValueError as exc:
            self.report.log(logging.ERROR, exc.args[0])
            # throw a fatal error, but assign version numbers to
            # see if more errors are caught
            major = 0
            minor = 0

        # check for previous uploads
        previous_uploads = ExperimentUpload.objects.filter(
            target=self.target,
            project=self.project,
        )
        if previous_uploads.exists():
            last_upload = previous_uploads.order_by("upload_version").last()

            data_version_validated, ver_val_msg = validate_data_version(
                major,
                minor,
                o_major=last_upload.data_version_major,
                o_minor=last_upload.data_version_minor,
            )
            upload_version_validated, upload_val_msg = validate_upload_version(
                self.version_number,
                previous_version=last_upload.upload_version,
            )
        else:
            data_version_validated, ver_val_msg = validate_data_version(major, minor)
            upload_version_validated, upload_val_msg = validate_upload_version(
                self.version_number,
                previous_version=0,  # no previous versions, exploit the +1 check
            )

        if not data_version_validated:
            self.report.log(logging.ERROR, ver_val_msg)

        if data_version_validated and ver_val_msg:
            self.report.log(logging.WARNING, ver_val_msg)

        if not upload_version_validated:
            self.report.log(logging.ERROR, upload_val_msg)

        if upload_version_validated and upload_val_msg:
            self.report.log(logging.WARNING, upload_val_msg)

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
        self.experiment_upload.data_version_major = major
        self.experiment_upload.data_version_minor = minor
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

        # just before actually processing objects, deprecate old
        # crystals
        # likely just a few objects at a time, if any
        for exp in Experiment.filter_manager.by_target(
            target=self.target,
        ).filter(
            code__in=self.excluded_crystals,
        ):
            exp.status = ExperimentStatusType.objects.get(status_code=4)
            exp.save()

        experiment_objects = self.process_experiment(
            yaml_data=crystals, prefix_tooltips=prefix_tooltips
        )
        # for val in experiment_objects.values():  # pylint: disable=no-member
        #     if val.new:
        #         val.instance.experiment_upload = self.experiment_upload
        #         val.instance.save()
        #         val.instance.refresh_from_db()

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

        # generate canon site objects and enumerate them starting from
        # where it was left off in the db. This is only needed for
        # tagging, but doing this here because don't want to pollute
        # tag module code with processed object dicts, trying to get
        # rid of them entirely
        canon_site_objects = self.process_canon_site(yaml_data=canon_sites)
        canon_sort_qs = (
            CanonSite.objects.filter(
                pk__in=[
                    k.instance.pk for k in canon_site_objects.values()
                ]  # pylint: disable=no-member
            )
            .annotate(
                # obvs=Count("canonsiteconf_set__siteobservation_set", default=0),
                obvs=Count("canonsiteconf__siteobservation", empty_result_set_value=0),
            )
            .order_by("-obvs", "name")
        )
        _canon_site_objects = {}
        for site in canon_sort_qs:
            key = f"{site.name}+{site.version}"
            _canon_site_objects[key] = canon_site_objects[
                key
            ]  # pylint: disable=no-member
        self._enumerate_objects(_canon_site_objects, "canon_site_num")
        for val in canon_site_objects.values():  # pylint: disable=no-member
            # instances modified, and will be modified down the line, refresh
            val.instance.refresh_from_db()

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
                        try:
                            exp_code_no = strip_exp_code(so.experiment.code)
                        except ValueError as exc:
                            self.report.log(logging.ERROR, exc.args[1])
                            # error, loading failed, use full code for demo
                            exp_code_no = so.experiment.code

                        code = f"{code_prefix}{exp_code_no}{next(suffix)}"

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
                        if code_qs.exists() and so.experiment.status.status_code == 0:
                            msg = (
                                f"short code {code} already exists for this target; "
                                + "specify a code_prefix to resolve this conflict"
                            )
                            self.report.log(logging.ERROR, msg)

                    so.code = code
                    so.save()

        for val in site_observation_objects.values():  # pylint: disable=no-member
            # instances modified, and will be modified down the line, refresh
            val.instance.refresh_from_db()

        # to be used in tagging, a necessity after the data exclusion (1674)
        so_qs = SiteObservation.filter_manager.by_target(self.target)

        # site_observations_versioned = {}
        # for val in site_observation_objects.values():  # pylint: disable=no-member
        #     site_observations_versioned[val.versioned_key] = val.instance

        # final remaining fk, attach reference site observation to canon_site_conf
        for val in canon_site_conf_objects.values():  # pylint: disable=no-member
            try:
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
            except KeyError as exc:
                # data may be missing. check the database
                try:
                    longcode = longcode_from_tag(val.index_data["reference_ligands"])
                    so = so_qs.get(longcode=longcode)
                    val.instance.ref_site_observation = so
                    val.instance.save()
                    msg = (
                        f"SiteObservation {val.index_data['reference_ligands']}"
                        + f" missing from {METADATA_FILE}"
                        + ", fetching from db",
                    )
                    logger.info(msg)
                except SiteObservation.DoesNotExist:
                    self.report.log(
                        logging.ERROR,
                        f"SiteObservation {val.index_data['reference_ligands']}"
                        + " missing from database",
                    )
            val.instance.refresh_from_db()

        # enumerate xtalform_sites. a bit trickier than others because
        # requires alphabetic enumeration starting from the letter of
        # the chain and following from there

        # sort the dictionary
        # fmt: off
        xtls_sort_qs = XtalformSite.objects.filter(
            pk__in=[k.instance.pk for k in xtalform_sites_objects.values() ], # pylint: disable=no-member
        ).annotate(
            obvs=Count("canon_site__canonsiteconf__siteobservation", empty_result_set_value=0),
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

        logger.debug("data read and processed, adding tags")
        tagger = TagManager(self.target)
        tagger.add_tags_to_canon_sites(
            canon_site_pks=[
                k.instance.pk
                for k in canon_site_objects.values()  # pylint: disable=no-member
            ]
        )
        tagger.add_tags_to_conformer_sites(
            canon_site_conf_pks=[
                k.instance.pk
                for k in canon_site_conf_objects.values()  # pylint: disable=no-member
            ]
        )
        tagger.add_tags_to_quatassemblies(
            quatassembly_pks=[
                k.instance.pk
                for k in quat_assembly_objects.values()  # pylint: disable=no-member
            ]
        )
        tagger.add_tags_to_xtalforms(
            xtalform_pks=[
                k.instance.pk
                for k in xtalform_objects.values()  # pylint: disable=no-member
            ]
        )
        tagger.add_tags_to_xtalformsites(
            xtalformsite_pks=[
                k.instance.pk
                for k in _xtalform_sites_objects.values()  # pylint: disable=no-member
            ]
        )
        # tag all new observations, so that the curator can find and
        # re-pose them
        datestr = timezone.now().date().strftime('%Y-%m-%d')
        tagger.tag_observations(
            f"{self.version_dir} {datestr}",
            "",
            category=TagCategory.objects.get(category="Other"),
            site_observations=[
                k.instance
                for k in site_observation_objects.values()  # pylint: disable=no-member
                if k.new
            ],
            clean_ids=False,
        )

        # see comment in method body if anything needs to be further
        # added after this method
        self._refresh_poses(site_observation_objects)
        self._generate_poses()

        # import compound identifier file, if present
        alias_file_path = (
            Path(upload_dir).joinpath("extra_files").joinpath(CUSTOM_IDENTIFIER_FILE)
        )
        if alias_file_path.exists():
            self.import_compound_identifiers(alias_file_path)

        # TODO: remove
        for val in site_observation_objects.values():  # pylint: disable=no-member
            if val.new:
                self._assign_observation_quality_status(
                    val.instance,
                    # val.index_data["auto_build_score"],
                )

        if soakdb_path:
            self.process_soakdb(db_file=str(soakdb_path))

        if self.version_number > 1 and self.target.computedset_set.exists():
            self.link_compounds_to_computedmolecules(site_observation_objects)

        self.mol_coords_to_db(site_observation_objects)

    def import_compound_identifiers(self, alias_file_path):
        try:
            df = pd.read_csv(alias_file_path)
        except UnicodeDecodeError:
            self.report.log(
                logging.ERROR,
                f"Error reading {CUSTOM_IDENTIFIER_FILE}, unexpected format",
            )
            return

        key_cols = ["xtal", "ligand_name"]
        extended_key_cols = key_cols + ["compound_code"]
        non_idf_cols = extended_key_cols + ["compound_code_update"]

        identifiers_from_file = [k for k in df.columns if k not in non_idf_cols]

        # mypy is doing it's thing again
        if not self.target.alias_order:  # type: ignore[attr-defined]
            self.target.alias_order = identifiers_from_file  # type: ignore[attr-defined]
            self.target.save()  # type: ignore[attr-defined]

        identifiers_from_file = set(identifiers_from_file)  # type: ignore[assignment]

        # I think this is a bad idea, but it was explicitly in the spec
        identifier_types = set(
            CompoundIdentifierType.objects.values_list("name", flat=True)
        )
        new_identifiers = identifiers_from_file.difference(identifier_types)  # type: ignore[attr-defined]
        for identifier in new_identifiers:
            CompoundIdentifierType(name=identifier).save()

        # you'd think I could supply the compounds processed, but I need a queryset..
        compounds = Compound.objects.annotate(
            exp_code=F("experiment__code"),
        ).filter(
            experiment__code__in=df["xtal"],
        )

        # validate cols, compound code should be unchanged
        for _, row in df[extended_key_cols].iterrows():
            exp_code, ligand_name, compound_code = row
            compound = compounds.get(exp_code=exp_code, ligand_name=ligand_name)
            if compound.compound_code != compound_code:
                self.report.log(
                    logging.ERROR,
                    (
                        f"{exp_code}, {ligand_name}: 'compound_code' not allowed to change."
                        + " use 'compound_code_update' column instead."
                    ),
                )

        # but if the correct column is supplied, then update
        if "compound_code_update" in df.columns:
            for _, row in df.loc[
                df["compound_code_update"].notna(), non_idf_cols
            ].iterrows():
                exp_code, ligand_name, _, compound_code_update = row
                compound = compounds.get(exp_code=exp_code, ligand_name=ligand_name)
                compound.compound_code = compound_code_update
                compound.save()

        identifiers = CompoundIdentifierType.objects.all()

        for idf in identifiers_from_file:
            identifier = identifiers.get(name=idf)

            for _, row in df.loc[df[idf].notna(), key_cols + [idf]].iterrows():
                exp_code, ligand_name, name = row
                compound = compounds.get(exp_code=exp_code, ligand_name=ligand_name)

                try:
                    compound_identifier = CompoundIdentifier(
                        type=identifier,
                        compound=compound,
                        name=name,
                    )
                    compound_identifier.save()

                    # set the preferred identifier to first non-empty
                    if not compound.current_identifier:
                        compound.current_identifier = compound_identifier
                        compound.save()

                except IntegrityError as exc:
                    # most probably a duplicate
                    self.report.log(
                        logging.ERROR,
                        exc.args[0],
                    )

                # memo to self: I tried using bulk_create here but it kept
                # failing with a rather cryptic error message. worth
                # revisiting when django has been upgraded

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

            if pose.main_site_observation.superseded:
                new_main = (
                    SiteObservation.filter_manager.by_target(
                        self.target,
                    )
                    .filter(
                        experiment=pose.main_site_observation.experiment,
                        cmpd=pose.main_site_observation.cmpd,
                        xtalform_site=pose.main_site_observation.xtalform_site,
                        canon_site_conf=pose.main_site_observation.canon_site_conf,
                        seq_id=pose.main_site_observation.seq_id,
                        chain_id=pose.main_site_observation.chain_id,
                    )
                    .order_by(
                        "-version",
                    )
                    .first()
                )

                pose.main_site_observation = new_main
                pose.save()

    def _refresh_poses(self, site_observation_objects):
        """Assign new main_observation if existing one has been superseded

        This is ran only on new site observation instances and
        *before* the pose generation, this way it skips the user
        modifications to poses.

        """

        for val in site_observation_objects.values():  # pylint: disable=no-member
            if val.new:
                logger.debug(
                    "processing poses for observation %s, %s, %s",
                    val.instance.pk,
                    val.instance.code,
                    val.instance.longcode,
                )
                # fmt: off
                qs = SiteObservation.filter_manager.by_target(
                    self.target,
                ).filter(
                    experiment=val.instance.experiment,
                    cmpd=val.instance.cmpd,
                    xtalform_site=val.instance.xtalform_site,
                    canon_site_conf=val.instance.canon_site_conf,
                    seq_id=val.instance.seq_id,
                    chain_id=val.instance.chain_id,
                    superseded=True,
                ).order_by(
                    "-version",
                )
                # fmt: on
                # older version(s) exist
                if qs.exists():
                    previous_main = qs.first()

                    # assign pose to new instance
                    val.instance.pose = previous_main.pose
                    val.instance.save()

                    # and then set the pose's main
                    previous_main.pose.main_site_observation = val.instance
                    previous_main.pose.save()

        # NB! this updates instances in the db but *not* in the
        # site_observation_objects dict. This means if another method
        # later operates on the site_observation instances inside the
        # dict and saves them, the changes made here will be lost. atm
        # the method is run at the end of the main processing method
        # and nothing after that saves the instances so that's fine,
        # but if something else needs to edit the observations,
        # refresh_from_db needs to be called

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

    def _assign_observation_quality_status(self, site_observation) -> None:
        status = QualityStatusType.objects.get(status="NONE")

        SiteObservationQualityStatus(
            site_observation=site_observation,
            status=status,
            user=None,
            auto_assigned=True,
            main_status=False,
            comment="Created on load",
        ).save()

    def link_compounds_to_computedmolecules(
        self, site_observation_objects: dict[str, MetadataObject]
    ) -> None:
        """Link incoming SiteObservations to existing ComputedMolecules.

        Spec (scraped from github (
        issue https://github.com/m2ms/fragalysis-frontend/issues/1591)).

        - on upload_1, do nothing
        - subsequent uploads, fFor every new molecule, check for a RHS
          design with the same chemical structure of the soaked
          compound (compare flattened inchikeys)
        - when match found, enumerate all possible LHS-RHS compound
          links, but annotate them with an RMSD
        - look only for ComputedMolecules within the target scope

        In practice, there's now a model
        SiteObservationComputedMolecule, effectively a m2m table
        between SiteObservation and ComputedMolecule that also
        captures an alignment RMSD value.
        """
        logger.debug('+linking observations to computed molecules')

        sdf_root = Path(settings.MEDIA_ROOT).joinpath(
            settings.COMPUTED_SET_MEDIA_DIRECTORY
        )

        # NB! see comment about filter_manager in managers.py for
        # compound only fetching LHS upload compounds. I believe here
        # this is the desired behaviour
        # compounds = Compound.filter_manager.by_target(self.target)

        # ComputedMolecules can come from two places:
        # - linked to a previously uploaded Compound
        # - linked to a previously uploaded ComputedSet
        # computed_molecules = ComputedMolecule.objects.filter(
        #     Q(computed_set__target=self.target) | Q(compound__in=compounds),
        # )
        computed_so = SiteObservation.objects.filter(
            xtalform_site__xtalform__in=Experiment.objects.filter(
                experiment_upload__target=self.target,
            ).values('xtalform'),
            experiment__isnull=True,
        )

        logger.debug('computed_siteobservations: %s', computed_so)
        for val in site_observation_objects.values():  # pylint: disable=no-member
            if not val.new:
                continue

            logger.debug('new observation: %s', val.instance)

            # NB! files have not been moved yet, need the tempdir location
            molpath = Path(settings.MEDIA_ROOT).joinpath(
                self.raw_data,
                *Path(val.instance.ligand_mol.name).parts[2:],
            )
            if not molpath.exists():
                continue

            logger.debug('molpath still going: %s', molpath)
            mol = Chem.MolFromMolFile(str(molpath))
            flattened_mol = copy.deepcopy(mol)
            Chem.RemoveStereochemistry(flattened_mol)
            flat_inchi = Chem.inchi.MolToInchiKey(flattened_mol)

            logger.debug('flat inchi: %s', flat_inchi)

            # the way the cset_loader is set up, the linked compound
            # is guaranteed to have a flattened inchi key. This is
            # explicitly used to .get() the compound instance and if
            # not found, new one is created. Which isn't really ideal,
            # just a missing inchi key may lead to duplicates. TODO
            # new issue and iron this out?
            logger.debug(
                'compmol set: %s',
                computed_so.filter(cmpd__inchi_key=flat_inchi),
            )
            for compmol in computed_so.filter(cmpd__inchi_key=flat_inchi):
                logger.debug('compmol: %s', compmol)

                # this can't be right. I need to compare 3D structures
                # but whatevs, going to be obsoleted
                # https://github.com/m2ms/fragalysis-frontend/issues/1748#issuecomment-3517491215

                # pr is it? did he mean the atom distance calc in cset_uplo?
                # cmol = Chem.MolFromMolBlock(compmol.sdf_info)
                # Chem.RemoveStereochemistry(cmol)

                logger.debug(
                    'cmol_path: %s', sdf_root.joinpath(str(compmol.virtual_ligand_mol))
                )
                cmol = Chem.MolFromMolFile(
                    sdf_root.joinpath(str(compmol.virtual_ligand_mol))
                )
                logger.debug('compmol_obj: %s', cmol)

                rmsd = None
                try:
                    rmsd = Chem.rdMolAlign.GetBestRMS(mol, cmol)
                    logger.debug('rmsd: %s', rmsd)
                except RuntimeError as exc:
                    # protection against rdkit internal errors
                    msg = (
                        f"Failed to find alignment between {compmol.virtual_molecule_name} "
                        + f'and {val.instance.code}'
                    )
                    # log an error, but don't stop processing
                    logger.error(msg)
                    logger.error(exc)

                # there is a unique constraint on this model, but only
                # new observations are being linked, so cannot clash
                # with any existing ones here
                SiteObservationComputedSiteObservation(
                    site_observation=val.instance,
                    computed_site_observation=compmol,
                    rmsd=rmsd,
                ).save()
                logger.debug('saved connection')

    def mol_coords_to_db(self, site_observation_objects):
        logger.debug('+mol_coords_to_db')
        for val in site_observation_objects.values():  # pylint: disable=no-member
            if val.new and val.index_data['mol']:
                logger.debug('coords for %s, %s', val.instance.pk, val.instance.code)
                mol = val.index_data['mol']
                conf = mol.GetConformer()

                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    atom = AtomCoordinates(
                        site_observation=val.instance,
                        coords=list(pos),
                        atom_number=atom.GetAtomicNum(),
                    )
                    atom.save()

    def exp_data_from_soakdb(self, row_data):
        # data structure to map db fields to functions that extract
        # data from soakdb. dictionary value is most of the time
        # simple parser that takes a field name as argument. If
        # argument is None, means it's more complex funciton and it
        # knows the fields it needs to operate on
        soakdb_field_resolvers = {
            "cchalf_high_res_shell": (
                self._soakdb_float,
                "DataProcessingCChalfHigh",
            ),
            "cchalf_overall": (
                self._soakdb_float,
                "DataProcessingCChalfOverall",
            ),
            "completeness_high_res_shell": (
                self._soakdb_float,
                "DataProcessingCompletenessHigh",
            ),
            "completeness_overall": (
                self._soakdb_float,
                "DataProcessingCompletenessOverall",
            ),
            "crystal_mounting_result": (
                self._soakdb_text,
                "MountingResult",
            ),
            "data_collection_date": (
                self._soakdb_datetime,
                "DataCollectionDate",
            ),
            "data_collection_outcome": (
                self._soakdb_data_collection_outcome,
                None,
            ),
            "dataset": (
                self._soakdb_text,
                "CrystalName",
            ),
            "date_model_last_updated": (
                self._soakdb_datetime,
                "LastUpdated",
            ),
            "date_status_updated": (
                self._soakdb_datetime,
                "RefinementOutcomeDate",
            ),
            "date_refined": (
                self._soakdb_datetime,
                "RefinementDate",
            ),
            "dimple_rfree": (
                self._soakdb_float,
                "DimpleRfree",
            ),
            "dimple_rwork": (
                self._soakdb_float,
                "DimpleRcryst",
            ),
            "experiment_comments": (
                self._soakdb_text,
                "SoakDBComments",
            ),
            "experiment_status": (
                self._soakdb_experiment_status,
                None,
            ),
            "experiment_type": (
                self._soakdb_experiment_type,
                None,
            ),
            "experiment_start_date": (
                self._soakdb_datetime,
                "SoakTimestamp",
            ),
            "final_compound_concentration_mm": (
                self._soakdb_float,
                "CompoundConcentration",
            ),
            "high_resolution": (
                self._soakdb_float,
                "DataProcessingResolutionHigh",
            ),
            "isig_i_overall": (
                self._soakdb_float,
                "DataProcessingIsigOverall",
            ),
            "isig_i_high_res_shell": (
                self._soakdb_float,
                "DataProcessingIsigHigh",
            ),
            "library": (
                self._soakdb_text,
                "LibraryName",
            ),
            "library_plate": (
                self._soakdb_text,
                "LibraryPlate",
            ),
            "ligand_confidence": (
                self._soakdb_ligand_confidence,
                None,
            ),
            "ligand_correlation_coefficient": (
                self._soakdb_text,
                "RefinementLigandCC",
            ),
            "model_last_updated_by": (
                self._soakdb_user,
                "LastUpdated_by",
            ),
            "modelled_smiles": (
                self._soakdb_text,
                "CompoundSMILES",
            ),
            "panddarun": (
                self._soakdb_text,
                "DimplePANDDApath",
            ),
            "pdb_code": (
                self._soakdb_text,
                "Deposition_PDB_ID",
            ),
            "processing_pipeline": (
                self._soakdb_text,
                "DataProcessingProgram",
            ),
            "refined_by": (
                self._soakdb_user,
                "RefinementRefiner",
            ),
            "refinement_comment": (
                self._soakdb_text,
                "RefinementComment",
            ),
            "refinement_rfree": (
                self._soakdb_float,
                "RefinementRfree",
            ),
            "refinement_rwork": (
                self._soakdb_float,
                "RefinementRcryst",
            ),
            "soakdb_entry": (
                self._soakdb_soakdb_entry,
                None,
            ),
            "soaking_time": (
                self._soakdb_duration,
                "SoakingTime",
            ),
            "source_well": (
                self._soakdb_text,
                "SourceWell",
            ),
            "space_group": (
                self._soakdb_space_group,
                None,
            ),
            "unit_cell_dimensions": (
                self._soakdb_numeric_array,
                "DataProcessingUnitCell",
            ),
            "refinement_resolution": (
                self._soakdb_float,
                "RefinementResolution",
            ),
        }
        exp_data = {}
        for db_field, (func, soakdb_field) in soakdb_field_resolvers.items():
            exp_data[db_field] = func(row_data, soakdb_field=soakdb_field)

        return exp_data

        # TODO: add other fields

    def process_soakdb(self, db_file: str) -> None:
        # fields to fetch from soakdb
        soakdb_fields = [
            "CompoundConcentration",
            "RefinementOutcome",
            "DataProcessingProgram",
            "DataProcessingCompletenessHigh",
            "SoakingTime",
            "RefinementRfree",
            "SoakTimestamp",
            "DataProcessingCChalfHigh",
            "RefinementRcryst",
            "LastUpdated_by",
            "LibraryPlate",
            "LastUpdated",
            "RefinementRefiner",
            "DataCollectionDate",
            "CrystalName",
            "DataProcessingIsigOverall",
            "DataProcessingIsigHigh",
            "LibraryName",
            "RefinementDate",
            "DimplePANDDApath",
            "Deposition_PDB_ID",
            "DataProcessingCompletenessOverall",
            "ID",
            "MountingResult",
            "DataProcessingCChalfOverall",
            "SourceWell",
            "RefinementSpaceGroup",
            "DataCollectionOutcome",
            "DimpleRcryst",
            "RefinementLigandCC",
            "RefinementOutcomeDate",
            "DimpleRfree",
            "DataProcessingResolutionHigh",
            "DataProcessingSpaceGroup",
            "DataProcessingUnitCell",
            "RefinementComment",
            "LabVisit",
            "CompoundSMILES",
            "SoakDBComments",
            "RefinementLigandConfidence",
            "RefinementResolution",
        ]

        query = f"SELECT {', '.join(soakdb_fields)} from mainTable"
        self.report.log(logging.INFO, "Processing soakdb")
        experiments: list[Experiment] = []
        qs = Experiment.filter_manager.by_target(self.target)
        with sqlite3.connect(db_file) as conn:
            conn.row_factory = sqlite3.Row  # This makes rows act like dicts
            cursor = conn.cursor()
            cursor.execute(query)
            rows = cursor.fetchall()
            logger.debug("Processing soakdb results")
            for row in rows:
                exp_data = self.exp_data_from_soakdb(dict(row))

                logger.debug("Extracted data for %s", exp_data["dataset"])
                # add some fields experiment needs to have
                exp_data["code"] = exp_data["dataset"]
                del exp_data["dataset"]

                exp_qs = qs.filter(code=exp_data["code"])
                if exp_qs.exists():
                    logger.debug("updating existing experiment: %s", exp_qs)
                    exp_qs.update(**exp_data)
                    exp = exp_qs.first()
                else:
                    logger.debug("Creating new experiment: %s", exp_data["code"])
                    exp_data["experiment_upload"] = self.experiment_upload
                    exp = Experiment(**exp_data)

                try:
                    exp.save()
                except IntegrityError as exc:
                    self.report.log(
                        logging.ERROR,
                        f"Failed to save experiment {exp_data['code']}: {str(exc)}",
                    )

                experiments.append(exp_data)

        self.report.log(
            logging.INFO, f"Processed {len(experiments)} entries from soakdb"
        )

    # functions extracting data from SoakDB. All have the same
    # signature, they take the dictionary of row header:value and
    # return a single value
    # CREATE DEFINER=`root`@`%` FUNCTION `soakDB`.`FUNCDataCollectionOutcome`(
    #     MountingResult VARCHAR(255),
    #     DataCollectionOutcome VARCHAR(255)
    # ) RETURNS varchar(255) CHARSET utf8mb4
    #     DETERMINISTIC
    # BEGIN

    # 	-- Determine data collection outcome based on the two input columns
    # 	IF MountingResult LIKE '%FAIL%' THEN
    # 		RETURN 'Failed - crystal did not survive soak';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%none%' THEN
    # 		RETURN 'Failed - autoprocessing failure';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%success%' THEN
    # 		RETURN 'Success - data collected';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%failed - low resolution%' THEN
    # 		RETURN 'Failed - low resolution';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%failed - centring failed%' THEN
    # 		RETURN 'Failed - centring failed';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%failed - processing%' THEN
    # 		RETURN 'Failed - processing';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%failed - no x-rays%' THEN
    # 		RETURN 'Failed - no X-rays';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%failed - loop broken%' THEN
    # 		RETURN 'Failed - loop broken';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%failed - unknown%' THEN
    # 		RETURN 'Failed - unknown';
    # 	ELSEIF DataCollectionOutcome IS NULL THEN
    # 		RETURN 'No data collected';
    #     ELSE
    #     	RETURN DataCollectionOutcome;
    #     END IF;

    # END;
    def _soakdb_data_collection_outcome(self, row_data, soakdb_field=None):
        del soakdb_field
        collection_outcome = row_data["DataCollectionOutcome"]
        mounting_result = row_data["MountingResult"]
        if mounting_result and mounting_result.find("FAIL") > -1:
            return "Failed - crystal did not survive soak"
        if collection_outcome:
            if collection_outcome.lower().find("nonw") > -1:
                return "Failed - autoprocessing failure"
            if collection_outcome.lower().find("success") > -1:
                return "Success - data collected"
            if collection_outcome.lower().find("failed - low resolution") > -1:
                return "Failed - low resolution"
            if collection_outcome.lower().find("failed - centring failed") > -1:
                return "Failed - centring failed"
            if collection_outcome.lower().find("failed - processing") > -1:
                return "Failed - processing"
            if collection_outcome.lower().find("failed - no x-rays") > -1:
                return "Failed - no X-rays"
            if collection_outcome.lower().find("failed - loop broken") > -1:
                return "Failed - loop broken"
            if collection_outcome.lower().find("failed - unknown") > -1:
                return "Failed - unknown"
            else:
                return collection_outcome
        else:
            return "No data collected"

    # CREATE DEFINER=`root`@`%` FUNCTION `soakDB`.`FUNCExperimentStatus`(
    #     DataCollectionOutcome VARCHAR(255),
    #     SoakTimestamp VARCHAR(20),
    #     RefinementOutcome VARCHAR(255)
    # ) RETURNS varchar(255) CHARSET utf8mb4
    #     DETERMINISTIC
    # BEGIN

    # 	-- Determine experiment status based on the three input columns
    # 	IF LOWER(DataCollectionOutcome) LIKE '%fail%' THEN
    # 		RETURN 'Data collection failed';
    # 	ELSEIF LOWER(DataCollectionOutcome) LIKE '%no data collected%' THEN
    # 		RETURN 'No data collected';
    # 	ELSEIF RefinementOutcome LIKE '%7 - Analysed & Rejected%' THEN
    # 		RETURN 'Structure analysed and rejected';
    # 	ELSEIF RefinementOutcome LIKE '%6 - Deposited%' THEN
    # 		RETURN 'Structure deposited';
    # 	ELSEIF RefinementOutcome LIKE '%5 - Deposition ready%' THEN
    # 		RETURN 'Deposition ready';
    # 	ELSEIF RefinementOutcome LIKE '%4 - CompChem ready%' THEN
    # 		RETURN 'CompChem ready';
    # 	ELSEIF RefinementOutcome LIKE '%3 - In Refinement%' THEN
    # 		RETURN 'Refinement in progress';
    # 	ELSEIF RefinementOutcome LIKE '%2 - PANDDA model%' THEN
    # 		RETURN 'Analysis in progress';
    # 	ELSEIF RefinementOutcome LIKE '%1 - Analysis Pending%' THEN
    # 		RETURN 'Analysis pending';
    # 	ELSEIF RefinementOutcome LIKE '%0 - All datasets%' THEN
    # 		RETURN 'Data collection failed';
    # 	ELSEIF RefinementOutcome IS NULL THEN
    # 		RETURN 'Data collection failed';
    # 	ELSEIF SoakTimestamp IS NULL THEN
    # 		RETURN 'Experiment pending';
    # 	ELSE
    # 		RETURN NULL;
    # 	END IF;

    # END;
    def _soakdb_experiment_status(self, row_data, soakdb_field=None):
        del soakdb_field
        collection_outcome = row_data["DataCollectionOutcome"]
        refinement_outcome = row_data["RefinementOutcome"]

        # this case isn't directly handled in trigger function, maybe
        # it's caught by final else
        if not collection_outcome:
            return None
        if not refinement_outcome:
            return "Data collection failed"
        if collection_outcome.lower().find("fail") > -1:
            return "Data collection failed"
        elif collection_outcome.lower().find("no data collected") > -1:
            return "No data collected"
        elif refinement_outcome.find("7 - Analysed & Rejecte") > -1:
            return "Structure analysed and rejected"
        elif refinement_outcome.find("6 - Deposited") > -1:
            return "Structure deposited"
        elif refinement_outcome.find("5 - Deposition ready") > -1:
            return "Structure deposited"
        elif refinement_outcome.find("4 - CompChem ready") > -1:
            return "CompChem ready"
        elif refinement_outcome.find("3 - In Refinement") > -1:
            return "Refinement in progress"
        elif refinement_outcome.find("2 - PANDDA model") > -1:
            return "Analysis in progress"
        elif refinement_outcome.find("1 - Analysis Pending") > -1:
            return "Analysis pending"
        elif not row_data["SoakTimestamp"]:
            return "Experiment pending"
        else:
            return None

    # CREATE DEFINER=`root`@`%` FUNCTION `soakDB`.`FUNCExperimentType`(SoakDBComments VARCHAR(255)) RETURNS varchar(255) CHARSET utf8mb4
    #     DETERMINISTIC
    # BEGIN

    # 	-- Get experiment type
    # 	IF LOWER(SoakDBComments) IN ('co-crystallisation','cocrystallisation','co-cryst','co-xtal','co-xtallisation','cocrystal') THEN
    # 		RETURN 'Co-crystallisation';
    # 	ELSE
    #         RETURN 'Soak';
    #     END IF;

    # END;
    @staticmethod
    def _soakdb_experiment_type(row_data, soakdb_field=None):
        del soakdb_field
        if row_data["SoakDBComments"]:
            if row_data["SoakDBComments"].lower() in (
                "co-crystallisation",
                "cocrystallisation",
                "co-cryst",
                "co-xtal",
                "co-xtallisation",
                "cocrystal",
            ):
                return 'Co-crystallisation'
            else:
                return "Soak"
        else:
            return None

    # CREATE DEFINER=`root`@`%` FUNCTION `soakDB`.`FUNCLigandConfidence`(
    #     RefinementLigandConfidence VARCHAR(255)
    # ) RETURNS varchar(255) CHARSET utf8mb4
    #     DETERMINISTIC
    # BEGIN

    # 	-- Determine ligand confidence
    # 	IF RefinementLigandConfidence LIKE '%0 - no ligand present%' THEN
    # 		RETURN 'No ligand modelled';
    # 	ELSEIF RefinementLigandConfidence LIKE '%1 - Low Confidence%' THEN
    # 		RETURN 'Low confidence';
    # 	ELSEIF RefinementLigandConfidence LIKE '%2 - Correct ligand, weak density%' THEN
    # 		RETURN 'Correct ligand, but weak density';
    # 	ELSEIF RefinementLigandConfidence LIKE '%3 - Clear density, unexpected ligand%' THEN
    # 		RETURN 'Unexpected ligand';
    # 	ELSEIF RefinementLigandConfidence LIKE '%4 - High Confidence%' THEN
    # 		RETURN 'High confidence';
    # 	ELSE
    # 		RETURN NULL;
    # 	END IF;

    # END;
    @staticmethod
    def _soakdb_ligand_confidence(row_data, soakdb_field=None):
        del soakdb_field
        confidence = row_data["RefinementLigandConfidence"]
        if confidence:
            if confidence.lower().find("0 - no ligand present") > -1:
                return "No ligand modelled"
            elif confidence.lower().find("1 - Low Confidence") > -1:
                return "Low confidence"
            elif confidence.lower().find("2 - Correct ligand, weak density") > -1:
                return "Correct ligand, but weak density"
            elif confidence.lower().find("3 - Clear density, unexpected ligand") > -1:
                return "Unexpected ligand"
            elif confidence.lower().find("4 - High Confidence") > -1:
                return "High confidence"
            else:
                return None
        else:
            return None

    @staticmethod
    def _soakdb_user(row_data, soakdb_field=None):
        soakdb_user = row_data[soakdb_field]

        # first check if any of the known unknowns
        if not soakdb_user or soakdb_user.lower().strip() == "none":
            return None
        elif soakdb_user.lower().strip() == "unknown":
            return get_user_model().objects.get(pk=settings.ANONYMOUS_USER)

        try:
            user = get_user_model().objects.get(username=soakdb_user.strip())
            return user
        except get_user_model().DoesNotExist:
            # create a placeholder
            placeholder = get_user_model()(
                username=soakdb_user.strip(),
                is_superuser=False,
                is_staff=False,
            )
            placeholder.save()
            return placeholder

    def _soakdb_soakdb_entry(self, row_data, soakdb_field=None):
        del soakdb_field
        return f"{row_data['ID']}_{row_data['LabVisit']}"

    # CREATE DEFINER=`root`@`%` FUNCTION `soakDB`.`FUNCSpacegroup`(
    #     RefinementSpaceGroup VARCHAR(255),
    #     DataProcessingSpaceGroup VARCHAR(255)
    # ) RETURNS varchar(255) CHARSET utf8mb4
    #     DETERMINISTIC
    # BEGIN

    # 	-- Check if RefinementSpaceGroup is NULL
    # 	IF RefinementSpaceGroup IS NULL THEN
    # 		RETURN REPLACE(IFNULL(DataProcessingSpaceGroup, ''), ' ', '');
    # 	ELSE
    # 		RETURN REPLACE(RefinementSpaceGroup, ' ', '');
    # 	END IF;

    # END;
    def _soakdb_space_group(self, row_data, soakdb_field=None):
        del soakdb_field
        if row_data["RefinementSpaceGroup"]:
            return row_data["RefinementSpaceGroup"].replace(" ", "")
        else:
            if row_data["DataProcessingSpaceGroup"]:
                return row_data["DataProcessingSpaceGroup"].replace(" ", "")
            else:
                # not directly handled by trigger code
                return None

    def _soakdb_float(self, row_data, soakdb_field=None):
        if row_data[soakdb_field] and row_data[soakdb_field] != "None":
            try:
                return float(row_data[soakdb_field])
            except ValueError:
                msg = (
                    f"Expected float on SoakDb ID:{row_data['ID']} {soakdb_field} "
                    + f"but received {row_data[soakdb_field]}"
                )
                self.report.log(logging.WARNING, msg)
                return None
        else:
            return None

    def _soakdb_numeric_array(self, row_data, soakdb_field=None):
        if row_data[soakdb_field] and row_data[soakdb_field] != "None":
            try:
                return [float(k) for k in row_data[soakdb_field].split()]
            except ValueError:
                msg = (
                    f"Expected numeric array on SoakDb ID:{row_data['ID']} {soakdb_field} "
                    + f"but received {row_data[soakdb_field]}"
                )
                self.report.log(logging.WARNING, msg)
                return None
        else:
            return None

    def _soakdb_datetime(self, row_data, soakdb_field=None):
        if row_data[soakdb_field] and row_data[soakdb_field] != "None":
            try:
                return parse(row_data[soakdb_field])
            except ParserError:
                # sometimes dates are given as:
                # 2020-12-02_09-50-12.03
                # cleanup:
                s_clean = row_data[soakdb_field].replace('_', ' ')
                s_clean = s_clean.replace('-', ':')
                try:
                    return parse(s_clean)
                except ParserError:
                    # still nothing
                    msg = (
                        f"Expected datetime on SoakDb ID:{row_data['ID']} {soakdb_field} "
                        + f"but received {row_data[soakdb_field]}"
                    )
                    self.report.log(logging.WARNING, msg)
                    return None
        else:
            return None

    def _soakdb_duration(self, row_data, soakdb_field=None):
        """Parses a string like '01:09:43' into a timedelta.

        NB! occasionally some incoming values may have 'AM'
        appended. Strip that.
        """
        logger.debug('duraton value: %s', row_data[soakdb_field])
        if row_data[soakdb_field] and row_data[soakdb_field] != "None":
            parts = row_data[soakdb_field].split(":")
            parts = [int(re.sub(r"\D", "", p)) for p in parts]

            # Support HH:MM:SS or MM:SS
            if len(parts) == 3:
                hours, minutes, seconds = parts
            elif len(parts) == 2:
                hours = 0
                minutes, seconds = parts
            else:
                msg = (
                    f"Expected timedelta on SoakDb ID:{row_data['ID']} {soakdb_field} "
                    + f"but received {row_data[soakdb_field]}"
                )
                self.report.log(logging.WARNING, msg)
                return None

            return timedelta(hours=hours, minutes=minutes, seconds=seconds)
        else:
            return None

    def _soakdb_text(self, row_data, soakdb_field=None):
        return row_data[soakdb_field]


def check_decompress_progress(process, archive_path, update, frequency=1.0):
    """Track decompression progress of tar+pigz.

    Following whichever process is actually reading the archive file.
    If anything goes wrong (proc not accessible, fd disappears, etc.),
    progress falls back to 0% without interfering with decompression.
    """
    try:
        total = os.path.getsize(archive_path)
    except OSError:
        total = 1  # avoid div-by-zero if file unreadable

    archive_real = os.path.realpath(archive_path)
    pid = process.pid
    fd_to_watch = None

    while process.poll() is None:
        pos = 0  # fallback default

        try:
            # if no fd yet, search /proc/<pid>/fd for one
            if fd_to_watch is None:
                try:
                    for fd in os.listdir(f"/proc/{pid}/fd"):
                        try:
                            link = os.readlink(f"/proc/{pid}/fd/{fd}")
                            if os.path.realpath(link) == archive_real:
                                fd_to_watch = fd
                                break
                        except OSError:
                            continue
                except FileNotFoundError:
                    pass  # /proc/<pid> might have disappeared

            # still nothing, maybe pigz is the one reading, not tar
            if fd_to_watch is None:
                try:
                    with open(f"/proc/{pid}/cmdline", "rb") as f:
                        cmdline = f.read().decode(errors="ignore")
                    if "tar" in cmdline and "-I" in cmdline:
                        # find child process (likely pigz)
                        for child in os.listdir("/proc"):
                            if not child.isdigit():
                                continue
                            try:
                                with open(  # pylint: disable=unspecified-encoding
                                    f"/proc/{child}/stat"
                                ) as f:
                                    parts = f.read().split()
                                    ppid = int(parts[3])
                                if ppid == pid:
                                    pid = int(child)  # switch to child process
                                    break
                            except Exception:
                                continue
                except Exception:
                    pass

            # got fd, read its position
            if fd_to_watch is not None:
                try:
                    with open(  # pylint: disable=unspecified-encoding
                        f"/proc/{pid}/fdinfo/{fd_to_watch}"
                    ) as f:
                        for line in f:
                            if line.startswith("pos:"):
                                pos = int(line.split()[1])
                                break
                except FileNotFoundError:
                    pos = 0  # fd vanished between checks
        except Exception:
            pos = 0  # catch-all, never crash

        # avoid division by 0 in progress computation
        progress = min(pos / total, 1.0) if total > 0 else 0.0
        try:
            update(progress)
        except Exception:
            # another catch-all, avoid crashing the decompression process
            pass

        time.sleep(frequency)


def load_target(
    data_bundle,
    proposal_ref=None,
    user_id=None,
    task=None,
):
    with TemporaryDirectory(dir=settings.MEDIA_ROOT) as tempdir:
        target_loader = TargetLoader(
            data_bundle, proposal_ref, tempdir, user_id=user_id, task=task
        )

        # Decompression can take some time, so we want to report progress
        bundle_filename = os.path.basename(data_bundle)
        target_loader.report.log(logging.INFO, f"Decompressing '{bundle_filename}'")

        try:
            msg = f"Extracting bundle: {data_bundle}"
            logger.info("%s%s", target_loader.report.task_id, msg)

            process = subprocess.Popen(
                [
                    "tar",
                    "-I",
                    "pigz",
                    "-xf",
                    target_loader.bundle_path,
                    "-C",
                    target_loader.raw_data,
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            check_decompress_progress(
                process,
                target_loader.bundle_path,
                lambda p: target_loader.report.log(
                    logging.INFO, f"Decompressing: {p:.1%}"
                ),
                frequency=5.0,
            )
            process.wait()

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
