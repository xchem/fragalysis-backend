"""Handler for successful Job execution. Basically initiates a number of
celery tasks to handle the Job result.
"""
import os
import shlex
from datetime import datetime

from celery.utils.log import get_task_logger
from django.http import HttpResponse

from viewer.models import JobRequest
from viewer.tasks import (
    erase_compound_set_job_material,
    process_compound_set,
    process_compound_set_job_file,
    validate_compound_set,
)

logger = get_task_logger(__name__)


def job_success_handler(
    jr: JobRequest,
    transition_time: datetime,
) -> HttpResponse:
    """Logic to start Squonk Job results retrieval. Called normally from the
    JobCallback endpoint on successful Job completion of a Job. The code here
    was extracted from the JobCallbackView for issue 1649 so that it could also
    be used from the JobRequestView (that can be triggered by a manual refresh).
    """
    assert jr
    assert transition_time

    assert jr.job_status == 'SUCCESS'

    code: str = str(jr.code)
    if jr.upload_status != 'PENDING':
        logger.warning(
            '- code=%s upload_status=%s. Ignoring, already uploading?',
            code,
            jr.upload_status,
        )
        return HttpResponse(status=204)

    # Change of status and SUCCESS
    # - mark the job upload as 'started'
    jr.upload_status = 'STARTED'
    jr.save()

    # We look for and upload the Job's output file.
    # For now there must be an '--outfile' in the job info's 'command'.
    # Here we have hard-coded the expectations because the logic to identify the
    # command's outputs is not fully understood.
    # The command is a string that we split and search.
    job_output = ''
    jr_job_info_msg = jr.squonk_job_info['msg']
    command = jr_job_info_msg.get('command')
    command_parts = shlex.split(command)
    outfile_index = 0
    while (
        outfile_index < len(command_parts)
        and command_parts[outfile_index] != '--outfile'
    ):
        outfile_index += 1
    # Found '--command'?
    if (
        command_parts[outfile_index] == '--outfile'
        and outfile_index < len(command_parts) - 1
    ):
        # Yes ... the filename is the next item in the list
        job_output = command_parts[outfile_index + 1]
    job_output_path = f'/{os.path.dirname(job_output)}'
    job_output_filename = os.path.basename(job_output)

    logger.info('code=%s job_output_path="%s"', code, job_output_path)

    # Initiate an upload (and removal) of files from Squonk.
    # Which requires the linking of several tasks.
    # We star the process with 'process_compound_set_job_file'
    # with the path and filename already discoverd...
    task_params = {
        'jr_id': jr.id,
        'transition_time': transition_time,
        'job_output_path': job_output_path,
        'job_output_filename': job_output_filename,
    }
    task_upload = (
        process_compound_set_job_file.s(task_params)
        | validate_compound_set.s()
        | process_compound_set.s()
        | erase_compound_set_job_material.s(job_request_id=jr.id)
    ).apply_async()

    logger.info(
        '- code=%s started process_job_file_upload(%s) task_upload=%s',
        code,
        jr.id,
        task_upload,
    )

    return HttpResponse(status=204)
