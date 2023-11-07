"""Functions (for the celery tasks) that transfer (download) files from
Fragalysis to Squonk.
"""
import urllib.parse
import os
import shutil

from django.conf import settings
from rest_framework import status
from squonk2.dm_api import DmApi

from celery.utils.log import get_task_logger
from viewer.models import (
    SiteObservation,
    ComputedMolecule,
    JobFileTransfer
)

logger = get_task_logger(__name__)


def process_file_transfer(auth_token, job_transfer_id):
    """Check/create the file transfer

    Args:
        auth_token from the request
        job_transfer_id
    """

    logger.info('+ process_file_transfer(job_transfer_id=%s)', job_transfer_id)

    job_transfer = JobFileTransfer.objects.get(id=job_transfer_id)
    logger.info('+ job_transfer.squonk_project=%s', job_transfer.squonk_project)

    # This to pick up NULL values from the changeover to using compounds.
    if not job_transfer.compounds:
        job_transfer.compounds = []

    num_proteins_to_transfer = len(job_transfer.proteins)
    num_compounds_to_transfer = len(job_transfer.compounds)
    num_to_transfer = num_proteins_to_transfer + num_compounds_to_transfer
    idx = 0
    logger.info('+ num_to_transfer=%s (%s + %s)',
                num_to_transfer, num_proteins_to_transfer, num_compounds_to_transfer)

    # The base directory for the source of the files we are transferring?
    # We expect files to include a path relative to TARGET_LOADER_MEDIA_DIRECTORY
    FILE_ROOT = os.path.join(settings.MEDIA_ROOT, settings.TARGET_LOADER_MEDIA_DIRECTORY)
    logger.info('+ FILE_ROOT=%s', FILE_ROOT)

    # Build the Squonk2 Project directory where files will be placed
    # e.g. "/fragalysis-files/hjyx".
    target = job_transfer.target
    squonk_directory = os.path.join('/', settings.SQUONK2_MEDIA_DIRECTORY, job_transfer.sub_path)
    logger.info('+ Destination squonk_directory=%s', squonk_directory)

    # All the files (proteins or compounds) are provided using relative
    # paths from the media directory. So we can join the tow lists
    # and treat them the same
    all_filename_refs = job_transfer.proteins + job_transfer.compounds
    if all_filename_refs:
        logger.info('+ Collecting files...')
        file_list = []
        for filename_ref in all_filename_refs:
            # We need to decode the file reference,
            # it is likely to be URL encoded.
            filename = urllib.parse.unquote(filename_ref)
            logger.info('+ Collecting filename=%s (target=%s)', filename, target)
            # File is expected to exist in the media directory
            file_path = os.path.join(FILE_ROOT, filename)
            if not os.path.isfile(file_path):
                msg = f'No such protein file ({file_path})'
                logger.error(msg)
                raise RuntimeError(msg)
            file_list.append(file_path)

        logger.info('+ Found %s files', len(file_list))
        logger.info('+ Calling DmApi.put_unmanaged_project_files() [proteins]...')
        result = DmApi.put_unmanaged_project_files(
            auth_token,
            project_id=job_transfer.squonk_project,
            project_files=file_list,
            project_path=squonk_directory,
            force=True,
        )
        logger.debug(result)

        if result.success:
            idx += 1
            job_transfer.transfer_progress = idx * 100 / num_to_transfer
            job_transfer.save()
            logger.info('+ Transferred files')
        else:
            msg = f'File Transfer Failed (msg={result.msg})'
            logger.error(msg)
            raise RuntimeError(msg)


def validate_file_transfer_files(request):
    """Check the request and return a list of proteins and/or computed molecule objects

    Args:
        request
    Returns
        error dict
        list of validated proteins (SiteObservation)
        list of validated computed molecules (ComputedMolecule)
    """
    error  = {}
    proteins = []
    compounds = []

    if request.data['proteins']:

        # Get first part of protein code
        proteins_list = [p.strip().split(":")[0]
                         for p in request.data['proteins'].split(',')]
        logger.info('+ Given proteins=%s', proteins_list)

        proteins = []
        for code_first_part in proteins_list:
            site_obvs = SiteObservation.objects.filter(code__contains=code_first_part).values()
            if site_obvs.exists():
                proteins.append(site_obvs.first())
            else:
                error['message'] = 'Please enter valid protein code for' \
                                   + ': {} '.format(code_first_part)
                error['status'] = status.HTTP_404_NOT_FOUND
                return error, proteins, compounds

        if len(proteins) == 0:
            error['message'] = 'API expects a list of comma-separated protein codes'
            error['status'] = status.HTTP_404_NOT_FOUND
            return error, proteins, compounds

    if request.data['compounds']:

        # Get compounds
        compounds_list = [c.strip() for c in request.data['compounds'].split(',')]
        logger.info('+ Given compounds=%s', compounds_list)

        compounds = []
        for compound in compounds_list:
            comp = ComputedMolecule.objects.filter(name=compound).values()
            if comp.exists():
                compounds.append(comp.first())
            else:
                error['message'] = 'Please enter valid compound name for' \
                                   + ': {} '.format(compound)
                error['status'] = status.HTTP_404_NOT_FOUND
                return error, proteins, compounds

        if len(compounds) == 0:
            error['message'] = 'API expects a list of comma-separated compound names'
            error['status'] = status.HTTP_404_NOT_FOUND
            return error, proteins, compounds

    if proteins or compounds:
        return error, proteins, compounds
    else:
        error['message'] = 'A valid set of protein codes and/or a list of valid' \
                           ' compound names must be provided'
        error['status'] = status.HTTP_404_NOT_FOUND
        return error, proteins, compounds
