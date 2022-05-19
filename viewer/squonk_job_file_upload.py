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
from dm_api.dm_api import DmApi

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

>  <original SMILES>  (1) 
O

>  <ref_url>  (1) 
{ref_url}

$$$$
"""

# Expected Job parameter file suffix
_JOB_PARAM_FILE_SUFFIX = '_param.json'

_DISABLE_DMAPI = True


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
    ref_url = settings.SQUONK_UI_URL
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


def process_compound_set_file(jr_id, transition_time):
    """Check the DM project for the expected file(s) and upload them.
    This code also applies the parameters to the uploaded compound set file.
    The full path to the uploaded/modified file is returned.

    Args:
        The user's authentication token
        The JobRequest record ID

    Returns the fully-qualified path to the SDF file (which may not exist as a file).
    """

    logger.info('Processing job compound file (%s)...', jr_id)

    job_request = JobRequest.objects.get(id=jr_id)

    # Do we need to create the upload path?
    # This is used for this 'job' and is removed when the upload is complete
    # successful or otherwise.
    upload_sub_dir = os.path.join('squonk_upload', str(job_request.code))
    upload_dir = create_media_sub_directory(upload_sub_dir)

    # There must be a 'variables.output' in the job's specification.
    job_squonk_job_spec = json.loads(job_request.squonk_job_spec)
    job_output_filename = job_squonk_job_spec.get('variables', {}).get('outfile')

    # The files we plan to upload to...
    tmp_sdf_filename = os.path.join(upload_dir, 'job.sdf')
    tmp_param_filename = os.path.join(upload_dir, 'job_params.json')

    if not job_output_filename:
        logger.warning("Squonk job spec had no 'outfile' (%d)", jr_id)
    else:

        # We expect an outfile (".sdf")
        # and an "{outfile}_params.json" file.

        if not _DISABLE_DMAPI:
            # Get the SDF file...
            result = DmApi.get_unmanaged_project_file(access_token=auth_token,
                                                      project_id=job_request.squonk_project,
                                                      project_file=job_output_filename,
                                                      local_file=tmp_sdf_filename,
                                                      timeout_s=120)
            logger.debug('DmApi.get_unmanaged_project_file(%s) result=%s',
                         job_output_filename, result)
            if result.success:
                # Get the parameter file file...
                job_output_param_filename = os.path\
                    .splitext(job_output_filename)[0] + _JOB_PARAM_FILE_SUFFIX
                result = DmApi.get_unmanaged_project_file(access_token=auth_token,
                                                          project_id=job_request.squonk_project,
                                                          project_file=job_output_param_filename,
                                                          local_file=tmp_param_filename,
                                                          timeout_s=20)
                logger.debug('DmApi.get_unmanaged_project_file(%s) result=%s',
                             job_output_param_filename, result)

    # If we have an SDF and parameter file
    # apply the blank molecule to the uploaded file
    # and then insert new parameters as we write to the actual SDF file.
    # i.e. we put a blank molecule in tmp.sdf and then apply parameters
    # as we re-write to {outfile}.
    sdf_filename = os.path.join(upload_dir, os.path.split(job_output_filename)[1])
    if os.path.isfile(tmp_sdf_filename) and os.path.isfile(tmp_param_filename):

        logger.info('Processing %s and %s', tmp_sdf_filename, tmp_param_filename)
        logger.info('Generating %s', sdf_filename)

        # Insert our 'blank molecule' into the uploaded file...
        _insert_sdf_blank_mol(job_request, transition_time, tmp_sdf_filename)
        with open(tmp_param_filename, 'r') as param_file:
            params = json.loads(param_file.read())
            add_prop_to_sdf(tmp_sdf_filename, sdf_filename, params)

    else:
        logger.warning('Not processing. Either %s or %s is missing',
                       tmp_sdf_filename, tmp_param_filename)

    logger.info('Done job compound file (%s)', jr_id)

    # Return the name of the file (which will exist if all's gone well).
    return sdf_filename
