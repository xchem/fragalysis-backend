"""
squonk_job_file_request Functions for creating squonk Jobs.

"""
from urllib.parse import urljoin
import os
import json
import logging
import datetime

from django.conf import settings
from dm_api.dm_api import DmApi

from viewer.models import ( Target,
                            Snapshot,
                            JobRequest,
                            JobFileTransfer )
from viewer.utils import get_https_host

logger = logging.getLogger(__name__)


def check_squonk_active(request):
    """Call a Squonk API to check that Squonk can be reached.
    """
    logger.info('+ Ping')
    auth_token = request.session['oidc_access_token']

    logger.info('oidc token')
    logger.info(auth_token)

    result = DmApi.ping(auth_token)
    logger.debug(result)

    return result.success


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
                result = DmApi.get_job(auth_token, job['id'])
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

    squonk_job_name = request.data['squonk_job_name']
    target_id = request.data['target']
    snapshot_id = request.data['snapshot']
    squonk_project = request.data['squonk_project']
    squonk_job_spec = request.data['squonk_job_spec']

    job_transfers = JobFileTransfer.objects.filter(snapshot=snapshot_id)
    if not job_transfers:
        raise ValueError('Files must be transferred before a job can be queued')

    job_transfer = JobFileTransfer.objects.filter(snapshot=snapshot_id).latest('id')
    if job_transfer.transfer_status != 'SUCCESS':
        raise ValueError('Job Transfer not complete')

    job_request = JobRequest()
    job_request.squonk_job_name = squonk_job_name
    job_request.user = request.user
    job_request.snapshot = Snapshot.objects.get(id=snapshot_id)
    job_request.target = Target.objects.get(id=target_id)
    job_request.squonk_project = squonk_project
    job_request.squonk_job_spec = squonk_job_spec

    # Saving creates the uuid for the callback
    job_request.save()
    callback_url = urljoin(get_https_host(request), os.path.join('api/job_callback',
                                                                 str(job_request.code)))

    # Ensure that the callback url ends with a slash so that the PUT works from Squonk
    callback_url = callback_url + '/'
    # Used for identifying the run, set to the username + date.
    job_name = job_request.user.username + '-' + datetime.date.today().strftime('%Y-%m-%d')

    logger.info('squonk_job_spec')
    logger.info(json.loads(squonk_job_spec))
    logger.info('callback url')
    logger.info(callback_url)
    logger.info('job_name')
    logger.info(job_name)

    result = DmApi.start_job_instance(auth_token,
                                      job_request.squonk_project,
                                      job_name,
                                      callback_url=callback_url,
                                      generate_callback_token=True,
                                      specification=json.loads(squonk_job_spec),
                                      timeout_s=8)
    logger.debug(result)

    if result.success:
        job_request.squonk_job_info = result
        job_request.squonk_url_ext = settings.SQUONK_INSTANCE_API + str(result.msg['instance_id'])
        job_request.job_start_datetime = datetime.datetime.utcnow()
        job_request.save()
        return job_request.id, job_request.squonk_url_ext

    job_request.delete()
    raise ValueError(result.msg)
