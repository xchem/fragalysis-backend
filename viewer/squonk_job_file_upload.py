"""Functions (for the celery tasks) that transfer (upload) files from
Squonk to Fragalysis.

These functions are used when a Job has been successful. It used the Job Request
to determine the output file(s) and then use the DM-API to get them. The code also
gets the expected generated parameters and adjusts the corresponding files.

The JobRequest record's 'squonk_job_spec' field will have a 'variables' section
that contains an 'outfile' declaration (e.g. "abc-molecules.sdf".).
"""
import json
import os
import shutil

from django.conf import settings
from rest_framework import status
from squonk2.dm_api import DmApi

from celery.utils.log import get_task_logger
from viewer.models import JobRequest, User
from viewer.utils import (
    SDF_VERSION,
    add_prop_to_sdf,
    create_media_sub_directory
)

logger = get_task_logger(__name__)

# A "Blank" molecule.
# Inserted at the top of SDF files pulled from Squonk2.
# This provides a 'header' that is then used
# to add the parameters found in the corresponding '_params.json' file.
#
# Use this string and a dictionary of values: -
#   _SDF_BLANK_MOL_TEMPLATE.format(**dictionary)
_SDF_BLANK_MOL_TEMPLATE = """{sdf_version}
     RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <submitter_name>  (1) 
{submitter_name}

>  <submitter_email>  (1) 
{submitter_email}

>  <submitter_institution>  (1) 
Diamond Light Source

>  <generation_date>  (1) 
{generation_date}

>  <method>  (1) 
{method}

>  <ref_url>  (1) 
{ref_url}

$$$$
"""

# Expected Job parameter file suffix.
# Jobs are acceptable candidates if they produce a file with this extension.
_JOB_PARAM_FILE_SUFFIX = '.meta.json'


def _insert_sdf_blank_mol(job_request, transition_time, sdf_filename):

    # Do nothing if the first line of the file matches the version we're about to set.
    blank_present = False
    with open(sdf_filename, 'r') as in_file:
        line = in_file.readline()
        if line and line.startswith(SDF_VERSION):
            blank_present = True
    if blank_present:
        return

    # Compound set reference URL.
    # What's the https-prefixed URL to the instance?
    # The record's URL is relative to the API.
    ref_url = settings.SQUONK2_UI_URL
    if ref_url.endswith('/'):
        ref_url += job_request.squonk_url_ext
    else:
        ref_url += '/' + job_request.squonk_url_ext

    # Compound set method.
    # This is restricted to 50 characters atm.
    # it's 49 characters if we set it to 'sq2-{instance UUID}'
    # which also happens to be the trailing part of the URL.
    instance_uuid = os.path.split(job_request.squonk_url_ext)[1]
    method = f'sq2-{instance_uuid}'

    # The transition time is an ISO8601 string,
    # with the date on the left of the 'T'

    variables = {'sdf_version': SDF_VERSION,
                 'submitter_name': job_request.user.username,
                 'submitter_email': job_request.user.email,
                 'generation_date': transition_time.split('T')[0],
                 'method': method,
                 'ref_url': ref_url}
    blank_mol = _SDF_BLANK_MOL_TEMPLATE.format(**variables)
    tmp_filename = sdf_filename + '.tmp'
    with open(tmp_filename, 'w') as tmp_file:
        tmp_file.write(blank_mol)
        with open(sdf_filename, 'r') as in_file:
            for line in in_file:
                tmp_file.write(line)
    os.remove(sdf_filename)
    os.rename(tmp_filename, sdf_filename)


def get_upload_sub_directory(job_request):
    """Returns the media-relative directory for the upload job data.
    """
    return os.path.join('squonk_upload', str(job_request.code))


def process_compound_set_file(jr_id,
                              transition_time,
                              job_output_path,
                              job_output_filename):
    """Check the DM project for the expected file(s) and upload them.
    This code also applies the parameters to the uploaded compound set file.
    The full path to the uploaded/modified file is returned.

    Args:
        The user's authentication token
        The JobRequest record ID
        The Job's output file path (prefixed with '/' and relative to the project root)
        The Job's filename (pathless, and an SD-File)

    Returns the fully-qualified path to the SDF file (sdf_filename).
    The file will not exist as a file if the upload processing fails.
    """

    logger.info('Processing job compound file (%s)...', jr_id)

    jr = JobRequest.objects.get(id=jr_id)

    # The callback token (required to make Squonk API calls from the callback)
    # and instance ID are in JobRequest.
    jr_job_info_msg = jr.squonk_job_info[1]
    token = jr_job_info_msg.get('callback_token')
    instance_id = jr_job_info_msg.get('instance_id')
    logger.info("Squonk API token=%s", token)
    logger.info("Squonk API instance_id=%s", instance_id)

    logger.info("Expecting Squonk path='%s'", job_output_path)

    # Do we need to create the upload path?
    # This is used for this 'job' and is removed when the upload is complete
    # successful or otherwise.
    upload_sub_dir = get_upload_sub_directory(jr)
    upload_dir = create_media_sub_directory(upload_sub_dir)

    # The temporary files we plan to upload to
    # (from the expected job_output_path and job_output_filename).
    tmp_sdf_filename = os.path.join(upload_dir, 'job.sdf')
    tmp_param_filename = os.path.join(upload_dir, 'job.meta.json')

    # The ultimate file we expect to create (after processing the temporary files).
    # Made unique with the JobRequest Squonk Instance ID.
    # This has to have the form '<upload_dir>/<prefix>_<cs_name>.sdf' and the
    # 'cs_name' must be unique as it will be used as the name of the CompoundSet
    # that will be created. We'll set this to '<path>/job_job-<job.code>.sdf'.
    #
    # The 'cs_name' must be no longer than 50 characters,
    # 'job-<uuid4>' should be 40 characters long.
    sdf_filename = os.path.join(upload_dir, f'job_job-{job.code}.sdf')

    # We expect an outfile (".sdf")
    # and an "{outfile}_params.json" file.
    got_all_files = False

    # Get the parameter file...
    #         --------------
    job_output_param_filename = os.path\
        .splitext(job_output_filename)[0] + _JOB_PARAM_FILE_SUFFIX
    logger.info("Expecting Squonk job_output_param_filename='%s'...",
                job_output_param_filename)
    result = DmApi.get_unmanaged_project_file_with_token(
        token=token,
        project_id=jr.squonk_project,
        project_path=job_output_path,
        project_file=job_output_param_filename,
        local_file=tmp_param_filename,
        timeout_s=20)
    if not result.success:
        logger.warning('Failed to get SDF parameter file (path=%s filename=%s)',
                       job_output_path, job_output_param_filename)
    else:
        # Got the parameter file so now get the SDF...
        #                               -----------
        logger.info("Expecting Squonk job_output_filename='%s'...",
                    job_output_filename)
        result = DmApi.get_unmanaged_project_file_with_token(
            token=token,
            project_id=jr.squonk_project,
            project_path=job_output_path,
            project_file=job_output_filename,
            local_file=tmp_sdf_filename,
            timeout_s=120)
        if not result.success:
            logger.warning('Failed to get SDF (path=%s filename=%s)',
                           job_output_path, job_output_filename)
        else:
            # Both files pulled back.
            got_all_files = True
            # Delete the callback token, which is no-longer needed.
            # Don't care if this fails - the token will expire automatically
            # after a period of time.
#            _ = DmApi.delete_instance_token(
#                instance_id=instance_id,
#                token=token)

    if not got_all_files:
        logger.warning('Not processing. Either %s or %s is missing',
                       tmp_sdf_filename, tmp_param_filename)
        return sdf_filename

    param_size = os.path.getsize(tmp_param_filename)
    sdf_size = os.path.getsize(tmp_sdf_filename)
    logger.info('Uploaded parameters to %s [%s bytes]', tmp_param_filename, param_size)
    logger.info('Uploaded SDF to %s [%s bytes]', tmp_sdf_filename, sdf_size)

    # Uploaded files cannot be empty.
    if param_size == 0 or sdf_size == 0:
        logger.warning('Not processing. Either the param file or SDF is empty')
        return sdf_filename

    # Got both files (non-empty) if we get gere...

    # Apply the blank molecule to the uploaded file
    # and then insert new parameters as we write to the actual SDF file.
    # i.e. we put a blank molecule in tmp.sdf and then apply parameters
    # as we re-write to {outfile}.

    logger.info('Generating %s...', sdf_filename)

    # Insert our 'blank molecule' into the uploaded file...
    _insert_sdf_blank_mol(jr, transition_time, tmp_sdf_filename)

    # Insert annotations...
    # The param file is a Squonk ".meta.json" file.
    # It should have an 'annotations' list with a 'fields' map in the list.
    #
    #   {'annotation_version':
    #           '0.0.1',
    #           'created': '2022-08-24T09:03:36.973340',
    #           'description': 'Fragmenstein combine',
    #           'fields': {'DDG': {'active': True,
    #                              'description': 'Delta deta G',
    #                              'required': True,
    #                              'type': 'number'},
    #                      'smiles': {'active': True,
    #                                 'description': 'Molecule SMILES',
    #                                 'required': False,
    #                                 'type': 'string'}},
    #
    # We take every field as the key and the description as the value.
    # If anything goes wrong we erase the SD file and return
    params = {}
    with open(tmp_param_filename, 'r') as param_file:
        meta = json.loads(param_file.read())
        if 'annotations' not in meta or len(meta['annotations']) == 0:
            logger.warning('Not processing. No annotations in %s',
                           tmp_param_filename)
            return sdf_filename
        for annotation in meta['annotations']:
            if 'fields' in annotation:
                for key, value in annotation['fields'].items():
                    if 'description' not in value:
                        logger.warning('Annotation field %s in %s has no "description"',
                                       key, tmp_param_filename)
                    else:
                        params[key] = value['description']

    if params:
        logger.info('Extracted %s from %s - adding to %s...',
                    params, tmp_param_filename, tmp_sdf_filename)
        add_prop_to_sdf(tmp_sdf_filename, sdf_filename, params)
        logger.info('Generated job compound file %s (%s)', sdf_filename, jr_id)
    else:
        logger.info('Found no suitable fields in %s', tmp_param_filename)

    # Return the name of the file
    # (which will exist if we pulled two files and all went well).
    return sdf_filename
