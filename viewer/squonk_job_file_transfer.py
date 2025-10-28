"""Functions (for the celery tasks) that transfer (download) files from
Fragalysis to Squonk.
"""
import os
import urllib.parse
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import unquote

from celery.utils.log import get_task_logger
from django.conf import settings
from rest_framework import status
from squonk2.dm_api import DmApi

from viewer.models import JobFileTransfer, SiteObservation

logger = get_task_logger(__name__)


def process_file_transfer(auth_token, job_transfer_id):
    """Check/create the file transfer

    Args:
        auth_token from the request
        job_transfer_id
    """

    logger.info('+ Processing file transfer (id=%s)', job_transfer_id)

    job_transfer = JobFileTransfer.objects.get(id=job_transfer_id)
    logger.info(
        '+ Transfer (id=%s) squonk_project=%s',
        job_transfer_id,
        job_transfer.squonk_project,
    )

    # This to pick up NULL values from the changeover to using compounds.
    if not job_transfer.compounds:
        job_transfer.compounds = []

    num_proteins_to_transfer = len(job_transfer.proteins)
    num_compounds_to_transfer = len(job_transfer.compounds)
    num_to_transfer = num_proteins_to_transfer + num_compounds_to_transfer
    logger.info(
        '+ Transfer (id=%s) num_to_transfer=%s (%s + %s)',
        job_transfer_id,
        num_to_transfer,
        num_proteins_to_transfer,
        num_compounds_to_transfer,
    )

    # Build the Squonk2 Project directory where files will be placed
    # e.g. "/fragalysis-files/hjyx".
    target = job_transfer.target
    squonk_directory = os.path.join(
        '/', settings.SQUONK2_MEDIA_DIRECTORY, job_transfer.sub_path
    )
    logger.info(
        '+ Transfer (id=%s) squonk_directory=%s', job_transfer_id, squonk_directory
    )
    if all_filename_refs := job_transfer.proteins + job_transfer.compounds:
        logger.info('+ Collecting files (id=%s)', job_transfer_id)
        file_list = []
        for filename_ref in all_filename_refs:
            # We need to decode the file reference,
            # it is likely to be URL encoded.
            filename = urllib.parse.unquote(filename_ref)
            logger.info(
                '+ Collecting %s (target=%s) (id=%s)', filename, target, job_transfer_id
            )
            # File is expected to exist in the media directory
            file_path = os.path.join(settings.MEDIA_ROOT, filename)
            if not os.path.isfile(file_path):
                msg = f'No such file {file_path} (id={job_transfer_id})'
                logger.error(msg)
                raise RuntimeError(msg)
            file_list.append(file_path)

        logger.info('+ Found %s files (id=%s)', len(file_list), job_transfer_id)
        logger.info(
            '+ Calling DmApi.put_unmanaged_project_files() (id=%s)', job_transfer_id
        )
        result = DmApi.put_unmanaged_project_files(
            auth_token,
            project_id=job_transfer.squonk_project,
            project_files=file_list,
            project_path=squonk_directory,
            force=True,
        )
        logger.debug(result)

        if result.success:
            job_transfer.transfer_progress = 100
            job_transfer.save()
            logger.info('+ Transferred files (id=%s)', job_transfer_id)
        else:
            msg = f'File Transfer Failed (id={job_transfer_id}) (msg={result.msg})'
            logger.error(msg)
            raise RuntimeError(msg)


def validate_file_transfer_files(
    request,
) -> Tuple[Optional[Dict[str, str]], Optional[List[Path]], Optional[List[Path]]]:
    """Check the request and return a list of proteins and/or computed molecule file
    path references (paths relative to the media directory).

    We're given a request that contains potentially empty strings that
    consist of comma-separated URL-encoded "proteins", and "compounds",
    and a "target" (target ID).

    Each protein and compound is a path and file to a file that is relative to the media
    directory. We just need to ensure that a SiteObservation exists for each
    (there should only be one) and it belongs to the given target.

    The user is already validated against the Target so here we check the given
    protein and compound references exist, and they belong to the Target.

    Args:
        request
    Returns
        error dictionary
        list of validated proteins
        list of validated computed molecules
    """
    assert 'target' in request.data
    assert 'proteins' in request.data
    assert 'compounds' in request.data

    logger.info('+ Validating file transfer files ()...')

    target_id = int(request.data['target'])

    protein_files: List[Path] = []
    compound_files: List[Path] = []

    if request.data['proteins']:
        # Get first part of protein code
        protein_paths_and_files = [
            unquote(p.strip()) for p in request.data['proteins'].split(',')
        ]
        for protein_path_and_file in protein_paths_and_files:
            if protein_path_and_file.endswith('_apo-desolv.pdb'):
                if not (
                    s_ob := SiteObservation.objects.filter(
                        apo_desolv_file=protein_path_and_file
                    ).first()
                ):
                    return tfr_validation_error(
                        f'Unknown Protein: "{protein_path_and_file}"',
                        status.HTTP_404_NOT_FOUND,
                    )

                s_ob_target_id = s_ob.experiment.experiment_upload.target.id
                if s_ob_target_id == target_id:
                    protein_files.append(Path(protein_path_and_file))
                else:
                    return tfr_validation_error(
                        f'Protein does not belong to Target: "{protein_path_and_file}"'
                        f' SiteObservation target={s_ob_target_id}'
                        f' Given target={target_id}',
                        status.HTTP_400_BAD_REQUEST,
                    )

        logger.info(
            "- Validated proteins (SiteObservations) [%d]",
            len(protein_files),
        )

    if request.data['compounds']:
        compound_paths_and_files = [
            unquote(p.strip()) for p in request.data['compounds'].split(',')
        ]
        for compound_path_and_file in compound_paths_and_files:
            if not (
                s_ob := SiteObservation.objects.filter(
                    ligand_mol=compound_path_and_file
                ).first()
            ):
                return tfr_validation_error(
                    f'Unknown Compound: "{compound_path_and_file}"',
                    status.HTTP_404_NOT_FOUND,
                )

            s_ob_target_id = s_ob.experiment.experiment_upload.target.id
            if s_ob_target_id == target_id:
                compound_files.append(Path(compound_path_and_file))
            else:
                return tfr_validation_error(
                    f'Compound does not belong to Target: "{compound_path_and_file}"'
                    f' SiteObservation target={s_ob_target_id}'
                    f' Given target={target_id}',
                    status.HTTP_400_BAD_REQUEST,
                )

        logger.info(
            "- Validated compounds (SiteObservations) [%d]",
            len(compound_files),
        )

    if not protein_files and not compound_files:
        return tfr_validation_error(
            'A valid set of protein codes and/or a list of valid compound names must be provided',
            status.HTTP_400_BAD_REQUEST,
        )

    logger.info(
        "- Validated file transfer files (%d, %d)",
        len(protein_files),
        len(compound_files),
    )
    return None, protein_files, compound_files


def tfr_validation_error(
    error: str, status_code: int
) -> Tuple[Dict[str, Any], None, None]:
    """Returns the error and HTTP status code as a tuple for a response."""
    return {'message': error, 'status': status_code}, None, None
