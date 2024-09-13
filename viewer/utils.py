import fnmatch
import itertools
import json
import logging
import os
import re
import shutil
import string
import tempfile
import uuid
from pathlib import Path
from typing import Dict, Generator, Optional
from urllib.parse import urlparse

import numpy as np
import pandas as pd
from django.conf import settings
from django.contrib.auth.models import User
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.core.mail import send_mail
from django.db import IntegrityError, transaction
from django.db.models import F
from django.http import JsonResponse
from rdkit import Chem

from scoring.models import SiteObservationGroup, SiteObvsSiteObservationGroup

from .models import (
    SiteObservation,
    SiteObservationTag,
    SiteObvsSiteObservationTag,
    Target,
)

logger = logging.getLogger(__name__)

# Set .sdf file format version
# Used at the start of every SDF file.
SDF_VERSION = 'ver_1.2'

SDF_RECORD_SEPARATOR = '$$$$\n'

# The root of all files constructed by 'dicttocsv'.
# The directory must not contain anything but dicttocsv-generated files.
# It certainly must not be the root of the media directory or any other directory in it.
# Introduced during 1247 security review.
CSV_TO_DICT_DOWNLOAD_ROOT = os.path.join(settings.MEDIA_ROOT, 'downloads', 'dicttocsv')


def is_url(url: Optional[str]) -> bool:
    try:
        result = urlparse(url)
        return all([result.scheme, result.netloc])
    except (ValueError, AttributeError):
        return False


def word_count(text: Optional[str]) -> int:
    """Returns an 'approximate' word count."""
    return len(text.split()) if text else 0


def create_squonk_job_request_url(instance_id):
    """Creates the Squonk Instance API url from an instance ID (UUID)."""
    return settings.SQUONK2_INSTANCE_API + str(instance_id)


def create_media_sub_directory(sub_directory):
    """Creates a directory (or directories) in the MEDIA directory,
    returning the full path.
    """
    directory = os.path.join(settings.MEDIA_ROOT, sub_directory)
    os.makedirs(directory, exist_ok=True)

    return directory


def delete_media_sub_directory(sub_directory):
    """Removes a media sub-directory."""
    assert sub_directory
    assert len(sub_directory)
    directory = os.path.normpath(os.path.join(settings.MEDIA_ROOT, sub_directory))
    if not os.path.isdir(directory):
        # No such directory!
        return
    if not directory.startswith(settings.MEDIA_ROOT) or os.path.samefile(
        directory, settings.MEDIA_ROOT
    ):
        # Danger!
        return

    shutil.rmtree(directory)


def add_prop_to_sdf(sdf_file_in, sdf_file_out, properties):
    """Returns an SDF file with the requested property
    (a dictionary of keys and values) added.
    Note that this assumes the file starts with a blank molecule.

    SDF parameters are of format:
        >  <TransFSScore>  (1)
        0.115601
    In this case we'd provide the following properties dictionary...
        {"TransFSScore": "0.115601"}
    """

    found_separator = False
    with open(sdf_file_out, 'a', encoding='utf-8') as sdf_out:
        with open(sdf_file_in, 'r', encoding='utf-8') as sdf_in:
            while True:
                if line := sdf_in.readline():
                    if not found_separator and line == SDF_RECORD_SEPARATOR:
                        # Found first separator
                        # dump the parameters now
                        found_separator = True
                        for name, value in properties.items():
                            sdf_out.write(f'>  <{name}>  (1)\n')
                            sdf_out.write(f'{value}\n')
                            sdf_out.write('\n')
                    sdf_out.write(line)
                else:
                    break


def add_props_to_sdf_molecule(
    *, sdf_file: str, properties: Dict[str, str], molecule: str
):
    """Given an input SDF, a dictionary of string properties and a molecule
    this function inserts the properties at the end of the molecule's record,
    just before the record separator. A temporary file is used that then replaces the
    input file.
    """
    # Strategy...
    # Search the file for the Molecule.
    # Then move to the end of record, and insert the properties.
    found_molecule: bool = False
    written_properties: bool = False
    with tempfile.NamedTemporaryFile(mode='a', encoding='utf-8') as temp:
        with open(sdf_file, 'r', encoding='utf-8') as sdf_in:
            while True:
                if line := sdf_in.readline():
                    if not found_molecule:
                        if line.strip() == molecule:
                            found_molecule = True
                    elif line == SDF_RECORD_SEPARATOR and not written_properties:
                        # Found end of molecule
                        # dump the parameters now
                        for name, value in properties.items():
                            temp.write(f'>  <{name}>\n')
                            temp.write(f'{value}\n')
                            temp.write('\n')
                        written_properties = True
                    # Write the original line
                    temp.write(line)
                else:
                    break
        # Flush the temporary file and replace the original file
        temp.flush()
        shutil.copy(temp.name, sdf_file)


def add_prop_to_mol(mol_field, mol_file_out, value):
    """Returns a mol_file with the requested property added.
    Note that the only property that seems to work with a .mol file if you write it is: "_Name".
    """
    rd_mol = Chem.MolFromMolBlock(mol_field)
    # Set the new property in the property mol object
    rd_mol.SetProp("_Name", value)
    Chem.MolToMolFile(rd_mol, mol_file_out)


# TODO: this method may be deprecated, not an issue with new uploads
def clean_filename(filepath):
    """Return the "clean" version of a Django filename without the '_abcdefg_' that is
    created when a file is overwritten.
    """
    file_split = os.path.splitext(os.path.basename(filepath))
    if fnmatch.fnmatch(
        file_split[0],
        '*_[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'
        + '[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]',
    ):
        cleaned_filename = file_split[0][:-8] + file_split[1]
    else:
        cleaned_filename = os.path.basename(filepath)
    return cleaned_filename


def get_https_host(request):
    """Common enabler code for returning https urls

    This is to map links to HTTPS to avoid Mixed Content warnings from Chrome browsers
    SECURE_PROXY_SSL_HEADER is referenced because it is used in redirecting URLs - if
    it is changed it may affect this code.
    Using relative links will probably also work, but This workaround allows both the
    'download structures' button and the DRF API call to work.
    Note that this link will not work on local
    """
    return settings.SECURE_PROXY_SSL_HEADER[1] + '://' + request.get_host()


def handle_uploaded_file(path: Path, f):
    with open(path, "wb+") as destination:
        for chunk in f.chunks(4096):
            destination.write(chunk)


def dump_curated_tags(filename: str) -> None:
    # fmt: off
    curated_tags = SiteObservationTag.objects.filter(
        user__isnull=False,
    ).annotate(
        ann_target_name=F('target__title'),
    )
    users = User.objects.filter(
        pk__in=curated_tags.values('user'),
    )
    siteobs_tag_group = SiteObvsSiteObservationTag.objects.filter(
        site_obvs_tag__in=curated_tags.values('pk'),
    ).annotate(
        ann_site_obvs_longcode=F('site_observation__longcode')
    )

    site_obvs_group = SiteObservationGroup.objects.filter(
        pk__in=curated_tags.values('mol_group'),
    ).annotate(
        ann_target_name=F('target__title'),
    )

    site_obvs_obvs_group = SiteObvsSiteObservationGroup.objects.filter(
        site_obvs_group__in=site_obvs_group.values('pk'),
    ).annotate(
        ann_site_obvs_longcode=F('site_observation__longcode')
    )
    # fmt: on

    result = {}
    for qs in (
        users,
        curated_tags,
        siteobs_tag_group,
        site_obvs_group,
        site_obvs_obvs_group,
    ):
        if qs.exists():
            jq = JsonResponse(list(qs.values()), safe=False)
            # have to pass through JsonResponse because that knows how
            # to parse django db field types
            data = json.loads(jq.content)
            name = qs[0]._meta.label  # pylint: disable=protected-access
            result[name] = data

    with open(filename, 'w', encoding='utf-8') as writer:
        writer.write(json.dumps(result, indent=4))


def restore_curated_tags(filename: str) -> None:
    with open(filename, 'r', encoding='utf-8') as reader:
        content = json.loads(reader.read())

    # models have to be saved in this order:
    # 1) User
    # 1) SiteObservationGroup  <- target
    # 2) SiteObservationTag    <- target, user
    # 3) SiteObvsSiteObservationGroup <- siteobvs
    # 3) SiteObvsSiteObservationTag <- siteobvs

    # takes a bit different approach with target and user - if user is
    # missing, restores the user and continues with tags, if target is
    # missing, skips the tag. This seems logical (at least at the time
    # writing this): if target hasn't been added obviously user
    # doesn't care about restoring the tags, but user might be
    # legitimately missing (hasn't logged in yet, and somebody else is
    # uploading the data)

    targets = Target.objects.all()
    site_observations = SiteObservation.objects.all()

    try:
        with transaction.atomic():
            new_mol_groups_by_old_pk = {}
            new_tags_by_old_pk = {}
            new_users_by_old_pk = {}

            user_data = content.get(
                User._meta.label,  # pylint: disable=protected-access
                [],
            )
            for data in user_data:
                pk = data.pop('id')
                try:
                    user = User.objects.get(username=data['username'])
                except User.DoesNotExist:
                    user = User(**data)
                    user.save()

                new_users_by_old_pk[pk] = user

            so_group_data = content.get(
                SiteObservationGroup._meta.label,  # pylint: disable=protected-access
                [],
            )
            for data in so_group_data:
                try:
                    target = targets.get(title=data['ann_target_name'])
                except Target.DoesNotExist:
                    logger.warning(
                        'Tried to restore SiteObservationGroup for target that does not exist: %s',
                        data['ann_target_name'],
                    )
                    continue

                data['target'] = target
                pk = data.pop('id')
                del data['ann_target_name']
                del data['target_id']
                sog = SiteObservationGroup(**data)
                sog.save()
                new_mol_groups_by_old_pk[pk] = sog

            so_tag_data = content.get(
                SiteObservationTag._meta.label,  # pylint: disable=protected-access
                [],
            )
            for data in so_tag_data:
                try:
                    target = targets.get(title=data['ann_target_name'])
                except Target.DoesNotExist:
                    logger.warning(
                        'Tried to restore SiteObservationTag for target that does not exist: %s',
                        data['ann_target_name'],
                    )
                    continue
                data['target'] = target
                pk = data.pop('id')
                del data['ann_target_name']
                del data['target_id']
                if data['mol_group_id']:
                    data['mol_group_id'] = new_mol_groups_by_old_pk[
                        data['mol_group_id']
                    ].pk
                data['user'] = new_users_by_old_pk[data['user_id']]
                del data['user_id']
                tag = SiteObservationTag(**data)
                try:
                    with transaction.atomic():
                        tag.save()
                except IntegrityError:
                    # this is an incredibly unlikely scenario where
                    # tag already exists - user must have, before
                    # restoring the tags, slightly edited an
                    # auto-generated tag. I can update the curated
                    # fields, but given they're both curated at this
                    # point, I choose to do nothing, skip the tag
                    logger.error(
                        'Curated tag %s already exists, skipping restore', data['tag']
                    )
                    continue

                new_tags_by_old_pk[pk] = tag

            so_so_group_data = content.get(
                SiteObvsSiteObservationGroup._meta.label,  # pylint: disable=protected-access
                [],
            )
            for data in so_so_group_data:
                try:
                    site_obvs = site_observations.get(
                        longcode=data['ann_site_obvs_longcode']
                    )
                except SiteObservation.DoesNotExist:
                    logger.warning(
                        'Tried to restore SiteObvsSiteObservationGroup for site_observation that does not exist: %s',
                        data['ann_site_obvs_longcode'],
                    )
                    continue
                site_obvs = site_observations.get(
                    longcode=data['ann_site_obvs_longcode']
                )
                data['site_observation'] = site_obvs
                del data['id']
                del data['ann_site_obvs_longcode']
                del data['site_observation_id']
                data['site_obvs_group'] = new_mol_groups_by_old_pk[
                    data['site_obvs_group_id']
                ]
                del data['site_obvs_group_id']
                SiteObvsSiteObservationGroup(**data).save()

            so_so_tag_data = content.get(
                SiteObvsSiteObservationTag._meta.label,  # pylint: disable=protected-access
                [],
            )
            for data in so_so_tag_data:
                try:
                    site_obvs = site_observations.get(
                        longcode=data['ann_site_obvs_longcode']
                    )
                except SiteObservation.DoesNotExist:
                    logger.warning(
                        'Tried to restore SiteObvsSiteObservationTag for site_observation that does not exist: %s',
                        data['ann_site_obvs_longcode'],
                    )
                    continue
                data['site_observation'] = site_obvs
                del data['id']
                del data['ann_site_obvs_longcode']
                del data['site_observation_id']
                data['site_obvs_tag'] = new_tags_by_old_pk.get(
                    data['site_obvs_tag_id'], None
                )
                del data['site_obvs_tag_id']
                if data['site_obvs_tag']:
                    # tag may be missing if not restored
                    SiteObvsSiteObservationTag(**data).save()

    except IntegrityError as exc:
        logger.error(exc)


def alphanumerator(
    start_from: str = "", drop_first: bool = True
) -> Generator[str, None, None]:
    """Return alphabetic generator (A, B .. AA, AB...) starting from a specified point.

    drop_first - as per workflow, usually it's given the last letter
    of previous sequence so the the next in the pipeline should be
    start_from + 1. drop_first = False indicates this is not necessary
    and start_from will be the first the iterator produces

    """

    # since product requries finite maximum return string length set
    # to 10 characters. that should be enough for fragalysis (and to
    # cause database issues)
    generator = (
        "".join(word)
        for word in itertools.chain.from_iterable(
            itertools.product(string.ascii_lowercase, repeat=i) for i in range(1, 11)
        )
    )

    # Drop values until the starting point is reached
    if start_from is not None and start_from != '':
        start_from = start_from.lower()
        generator = itertools.dropwhile(lambda x: x != start_from, generator)  # type: ignore[assignment]
        if drop_first:
            # drop one more, then it starts from after the start from as it should
            _ = next(generator)

    return generator


def sanitize_boolean_column(column):
    """
    Sanitize a DataFrame column to boolean values.
    Handles various representations of booleans such as:
    - Python booleans
    - String representations ('True', 'False', 'true', 'false', 'yes', 'no')
    - Numeric representations (1, 0)
    """
    # don't understans why it considers int and str numbers to be the same
    true_values = {  # pylint: disable=duplicate-value
        'true',
        '1',
        't',
        'y',
        'yes',
        'True',
        'TRUE',
        True,
        1,
    }
    false_values = {  # pylint: disable=duplicate-value
        'false',
        '0',
        'f',
        'n',
        'no',
        'False',
        'FALSE',
        False,
        0,
        None,
        '',
        np.nan,
    }

    def convert_to_boolean(value):
        if value in true_values:
            return True
        elif value in false_values:
            return False
        else:
            raise ValueError(f'Unexpected boolean value: {value}')

    return column.apply(convert_to_boolean)


def save_tmp_file(myfile):
    """Save file in temporary location for validation/upload processing"""

    name = myfile.name
    path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
    return str(os.path.join(settings.MEDIA_ROOT, path))


def create_csv_from_dict(input_dict, title=None, filename=None):
    """Write a CSV file containing data from an input dictionary and return a full
    to the file (in the media directory).
    """
    if not filename:
        filename = 'download'

    unique_dir = str(uuid.uuid4())
    download_path = os.path.join(CSV_TO_DICT_DOWNLOAD_ROOT, unique_dir)
    os.makedirs(download_path, exist_ok=True)

    download_file = os.path.join(download_path, filename)

    # Remove file if it already exists
    if os.path.isfile(download_file):
        os.remove(download_file)

    with open(download_file, "w", newline='', encoding='utf-8') as csvfile:
        if title:
            csvfile.write(title)
            csvfile.write("\n")

    df = pd.DataFrame.from_dict(input_dict)
    df.to_csv(download_file, mode='a', header=True, index=False)

    return download_file


def email_task_completion(
    contact_email, message_type, target_name, target_path=None, task_id=None
):
    """Notify user of upload completion"""

    logger.info('+ email_notify_task_completion: ' + message_type + ' ' + target_name)
    email_from = settings.EMAIL_HOST_USER

    if contact_email == '' or not email_from:
        # Only send email if configured.
        return

    if message_type == 'upload-success':
        subject = 'Fragalysis: Target: ' + target_name + ' Uploaded'
        message = (
            'The upload of your target data is complete. Your target is available at: '
            + str(target_path)
        )
    elif message_type == 'validate-success':
        subject = 'Fragalysis: Target: ' + target_name + ' Validation'
        message = (
            'Your data was validated. It can now be uploaded using the upload option.'
        )
    else:
        # Validation failure
        subject = 'Fragalysis: Target: ' + target_name + ' Validation/Upload Failed'
        message = (
            'The validation/upload of your target data did not complete successfully. '
            'Please navigate the following link to check the errors: validate_task/'
            + str(task_id)
        )

    recipient_list = [
        contact_email,
    ]
    logger.info('+ email_notify_task_completion email_from: %s', email_from)
    logger.info('+ email_notify_task_completion subject: %s', subject)
    logger.info('+ email_notify_task_completion message: %s', message)
    logger.info('+ email_notify_task_completion contact_email: %s', contact_email)

    # Send email - this should not prevent returning to the screen in the case of error.
    send_mail(subject, message, email_from, recipient_list, fail_silently=True)
    logger.info('- email_notify_task_completion')
    return


def sanitize_directory_name(name: str, path: Path | None = None) -> str:
    """
    Sanitize a string to ensure it only contains characters allowed in UNIX directory names.

    Parameters:
    name: The input string to sanitize.
    path (optional): the parent directory where the directory would reside, to check if unique

    Returns:
    str: A sanitized string with only allowed characters.
    """
    # Define allowed characters regex
    allowed_chars = re.compile(r'[^a-zA-Z0-9._-]')

    # Replace disallowed characters with an underscore
    sanitized_name = allowed_chars.sub('_', name.strip())

    # Replace multiple underscores with a single underscore
    sanitized_name = re.sub(r'__+', '_', sanitized_name)
    logger.debug('sanitized name: %s', sanitized_name)
    if path:
        target_dirs = [d.name for d in list(path.glob("*")) if d.is_dir()]
        logger.debug('target dirs: %s', target_dirs)
        new_name = sanitized_name
        suf = 1
        while new_name in target_dirs:
            suf = suf + 1
            new_name = f'{sanitized_name}_{suf}'
            logger.debug('looping suffix: %s', new_name)

        sanitized_name = new_name

    return sanitized_name


def clean_object_id(name: str) -> str:
    """Replace '/' and '+' with '/' in XCA object identifiers"""
    splits = name.split('-x')
    if len(splits) > 1:
        return f"{splits[0]}-x{splits[1].replace('+', '/').replace('_', '/')}"
    else:
        return name.replace('+', '/').replace('_', '/')
