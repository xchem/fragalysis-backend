"""
squonk_job_file_request Functions for creating squonk Jobs.

"""
from urllib.parse import urljoin
import os
import json
import logging
import datetime

import shortuuid

from squonk2.dm_api import DmApi

from viewer.models import ( Target,
                            Snapshot,
                            JobRequest,
                            JobFileTransfer )
from viewer.utils import create_squonk_job_request_url, get_https_host
from viewer.squonk2_agent import CommonParams, Squonk2Agent, get_squonk2_agent

logger = logging.getLogger(__name__)

_SQ2A: Squonk2Agent = get_squonk2_agent()

def check_squonk_active(request):
    """Call the Squonk2 Agent to check that Squonk2 can be reached.
    """
    logger.info('+ Squonk2Agent.ping()')

    ping_rv = _SQ2A.ping()
    logger.debug(ping_rv)

    return ping_rv.success


def get_squonk_job_config(request,
                          job_collection=None,
                          job_name=None,
                          job_version=None):
    """get squonk job configuration details from squonk.
    1. Get all available jobs for the user.
    2. Filter the job for the specific job name, collection and version if provided

    Returns:
        DICT
        either the list of available jobs
        or details for the requested job.
    """

    logger.info('+ get_squonk_job_config')
    auth_token = request.session['oidc_access_token']
    logger.debug(auth_token)

    result = DmApi.get_available_jobs(auth_token)
    logger.debug(result)

    if result.success is True:
        available_jobs = result.msg
    else:
        return result.msg

    if not job_name:
        # No name provided - return all
        return available_jobs
    else:
        # Got a job name (and collection and version)
        for job in available_jobs['jobs']:
            if job['job'] == job_name\
                    and job['collection'] == job_collection \
                    and job['version'] == job_version:
                result = DmApi.get_job(auth_token, job_id=job['id'])
                # This either returns the definition or the squonk message.
                return result.msg
        return {'Job not found'
                f' (collection={job_collection} name={job_name} version={job_version}'}


def create_squonk_job(request):
    """Check set up and create Squonk job instance.
    1. Check files have been successfully transferred for the snapshot.
    2. Queue job

    Return:
        the job id
        a URL allowing the front end to link to the running job instance
    """

    logger.info('+ create_squonk_job')
    auth_token = request.session['oidc_access_token']
    logger.debug(auth_token)

    access_id = request.data['access']  # The access ID/legacy Project record title
    target_id = request.data['target']
    snapshot_id = request.data['snapshot']
    session_project_id = request.data['session_project'] 
    squonk_job_name = request.data['squonk_job_name']
    squonk_job_spec = request.data['squonk_job_spec']

    logger.info('+ access_id=%s', access_id)
    logger.info('+ target_id=%s', target_id)
    logger.info('+ snapshot_id=%s', snapshot_id)
    logger.info('+ session_project_id=%s', session_project_id)
    logger.info('+ squonk_job_name=%s', squonk_job_name)
    logger.info('+ squonk_job_spec=%s', squonk_job_spec)

    job_transfers = JobFileTransfer.objects.filter(snapshot=snapshot_id)
    if not job_transfers:
        logger.warning('No JobFileTransfer object for snapshot %s', snapshot_id)
        raise ValueError('No JobFileTransfer object for snapshot %s.'
                         ' Files must be transferred before a job can be requested.',
                         snapshot_id)

    job_transfer = JobFileTransfer.objects.filter(snapshot=snapshot_id).latest('id')
    if job_transfer.transfer_status != 'SUCCESS':
        logger.warning('JobFileTransfer status is not SUCCESS (is %s)',
                       job_transfer.transfer_status)
        raise ValueError('Job Transfer not complete')

    logger.info('+ Calling ensure_project() to get the Squonk2 Project...')

    # This requires a Squonk2 Project (created by the Squonk2Agent).
    # It may be an existing project, or it might be a new project.
    user = request.user
    common_params = CommonParams(user_id=user.id,
                                 access_id=access_id,
                                 target_id=target_id,
                                 session_id=session_project_id)
    sq2_rv = _SQ2A.ensure_project(common_params)
    if not sq2_rv.success:
        msg = f'JobTransfer failed to get/create a Squonk2 Project' \
                f' for User "{user.username}", Access ID {access_id},' \
                f' Target ID {target_id}, and SessionProject ID {session_project_id}.' \
                f' Got "{sq2_rv.msg}".' \
                ' Cannot continue'
        logger.error(msg)
        raise ValueError(msg)

    # The Squonk2Project record is in the response msg
    squonk2_project_uuid = sq2_rv.msg.uuid
    squonk2_unit_name = sq2_rv.msg.unit.name
    squonk2_unit_uuid = sq2_rv.msg.unit.uuid
    logger.info('+ ensure_project() returned Project uuid=%s (unit="%s" unit_uuid=%s)',
                squonk2_project_uuid, squonk2_unit_name, squonk2_unit_uuid)

    job_request = JobRequest()
    job_request.squonk_job_name = squonk_job_name
    job_request.user = request.user
    job_request.snapshot = Snapshot.objects.get(id=snapshot_id)
    job_request.target = Target.objects.get(id=target_id)
    # We should use a foreign key,
    # but to avoid migration issues with the existing code
    # we continue to use the project UUID string field.
    job_request.squonk_project = squonk2_project_uuid
    job_request.squonk_job_spec = squonk_job_spec

    # Saving creates the uuid for the callback
    job_request.save()
    callback_url = urljoin(get_https_host(request), os.path.join('api/job_callback',
                                                                 str(job_request.code)))

    # Ensure that the callback url ends with a slash so that the PUT works from Squonk
    callback_url = callback_url + '/'
    # Used for identifying the run, set to the username + date.
    job_name = job_request.user.username + '-' + datetime.date.today().strftime('%Y-%m-%d')

    # Create a callback token
    # that we can use on the job from our callback context.
    # It is required to be a shortuuid of 22 characters using the default character set
    callback_token = shortuuid.uuid()

    logger.info('+ job_name=%s', job_name)
    logger.info('+ callback_url=%s', callback_url)
    logger.info('+ callback_token=%s', callback_token)
    
    # Dry-run the Job execution (this provides us with the 'command', which is 
    # placed in the JobRecord's squonk_job_info).
    logger.info('+ Calling DmApi.dry_run_job_instance(%s)', job_name)
    result = DmApi.dry_run_job_instance(auth_token,
                                        project_id=job_request.squonk_project,
                                        name=job_name,
                                        callback_url=callback_url,
                                        callback_token=callback_token,
                                        specification=json.loads(squonk_job_spec))
    logger.debug(result)

    if not result.success:
        logger.warning('+ dry_run_job_instance(%s) result=%s', job_name, result)
        logger.error('+ FAILED (configuration problem) (%s)', job_name)
        job_request.delete()
        raise ValueError(result.msg)

    # We can now commit the JobRequest record so that it's ready
    # for use by any callbacks. The 'result' will contain the callback token
    # and the Job's decoded command (that will be interrogated when the Job s complete)
    job_request.squonk_job_info = result
    job_request.job_start_datetime = datetime.datetime.utcnow()
    job_request.save()

    # Now start the job 'for real'...
    logger.info('+ Calling DmApi.start_job_instance(%s)', job_name)
    result = DmApi.start_job_instance(auth_token,
                                      project_id=job_request.squonk_project,
                                      name=job_name,
                                      callback_url=callback_url,
                                      callback_token=callback_token,
                                      specification=json.loads(squonk_job_spec),
                                      timeout_s=8)
    logger.debug(result)

    if not result.success:
        logger.warning('+ start_job_instance(%s) result=%s', job_name, result)
        logger.error('+ FAILED (job probably did not start) (%s)', job_name)
        job_request.delete()
        raise ValueError(result.msg)

    job_instance_id = result.msg['instance_id']
    job_task_id = result.msg['task_id']
    logger.info('+ SUCCESS. Job "%s" started (job_instance_id=%s job_task_id=%s)',
                job_name, job_instance_id, job_task_id)

    # Manufacture the Squonk URL (actually set in the callback)
    # so the front-end can use it immediately (we cannot set the JobRequest
    # `squonk_url_ext` here as it introduces a race-condition with the callback logic).
    squonk_url_ext = create_squonk_job_request_url(job_instance_id)

    return job_request.id, squonk_url_ext
