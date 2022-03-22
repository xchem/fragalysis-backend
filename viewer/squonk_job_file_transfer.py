"""
squonk_job_file_transfer Functions for the celery tasks for transferring
files from Fragalysis to Squonk.

"""

from .squonk import put_file


def check_file_transfer(request,
                        target,
                        snapshot,
                        proteins):
    """Check/create the file transfer for list of proteins

    Args:
        request
        target
        snapshot
        proteins

    Returns:
        Boolean:
    """

    # Settings for if validate option selected
    # Start celery task

    return True
