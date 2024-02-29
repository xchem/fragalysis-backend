import logging

from api.security import ISpyBSafeStaticFiles, ISpyBSafeStaticFiles2
from viewer.models import SiteObservation, Target

logger = logging.getLogger(__name__)

# def prot_download(request, file_path):
#     """
#     Download a protein by nginx redirect
#     :param request: the initial request
#     :param file_path: the file path we're getting from the static
#     :return: the response (a redirect to nginx internal)
#     """
#     logger.info("+ Received file_download file path: %s", file_path)
#     ispy_b_static = ISpyBSafeStaticFiles2()
#     ispy_b_static.model = SiteObservation
#     ispy_b_static.request = request
#     # ispy_b_static.permission_string = "target_id__project_id"
#     # ispy_b_static.permission_string = "target__project_id"
#     ispy_b_static.permission_string = "experiment__experiment_upload__target__project_id"
#     # ispy_b_static.field_name = "pdb_info"
#     ispy_b_static.field_name = "apo_file"
#     ispy_b_static.content_type = "application/x-pilot"
#     ispy_b_static.prefix = "/pdbs/"
#     # ispy_b_static.prefix = "/target_loader_data/"
#     ispy_b_static.input_string = file_path
#     return ispy_b_static.get_response()


def file_download(request, file_path):
    """
    Download a protein by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    logger.info("+ Received file_download file path: %s", file_path)
    ispy_b_static = ISpyBSafeStaticFiles2()
    # ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = SiteObservation
    ispy_b_static.request = request
    # ispy_b_static.permission_string = "target_id__project_id"
    # the following 2 aren't used atm
    ispy_b_static.permission_string = (
        "experiment__experiment_upload__target__project_id"
    )
    # ispy_b_static.field_name = "pdb_info"
    ispy_b_static.field_name = "apo_file"
    ispy_b_static.content_type = "application/x-pilot"
    # ispy_b_static.prefix = "target_loader_data/48225dbf-204a-48e1-8ae7-f1632f4dba89/Mpro-v2/Mpro/upload_2/aligned_files/Mpro_Nterm-x0029/"
    # ispy_b_static.prefix = "target_loader_data"
    # ispy_b_static.prefix = "/target_loader_data/"
    ispy_b_static.prefix = "/pdbs/"
    ispy_b_static.input_string = file_path
    return ispy_b_static.get_response()


def tld_download(request, file_path):
    """
    Download a protein by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    logger.info("+ Received tld_download file path: %s", file_path)
    ispy_b_static = ISpyBSafeStaticFiles2()
    # ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = SiteObservation
    ispy_b_static.request = request
    # ispy_b_static.permission_string = "target_id__project_id"
    # the following 2 aren't used atm
    ispy_b_static.permission_string = (
        "experiment__experiment_upload__target__project_id"
    )
    ispy_b_static.field_name = "apo_file"
    ispy_b_static.content_type = "application/x-pilot"
    ispy_b_static.prefix = "/target_loader_data/"
    ispy_b_static.input_string = file_path
    return ispy_b_static.get_response()


def cspdb_download(request, file_path):
    """
    Download a protein by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    logger.info("+ Received cspdb_download file path: %s", file_path)
    ispy_b_static = ISpyBSafeStaticFiles2()
    ispy_b_static.model = SiteObservation
    ispy_b_static.request = request
    # the following 2 aren't used atm
    ispy_b_static.permission_string = (
        "experiment__experiment_upload__target__project_id"
    )
    ispy_b_static.field_name = "apo_file"
    ispy_b_static.content_type = "application/x-pilot"
    ispy_b_static.prefix = "/computed_set_data/"
    ispy_b_static.input_string = file_path
    return ispy_b_static.get_response()


def bound_download(request, file_path):
    """
    Download a protein by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = SiteObservation
    ispy_b_static.request = request
    ispy_b_static.permission_string = "target_id__project_id"
    ispy_b_static.field_name = "bound_info"
    ispy_b_static.content_type = "application/x-pilot"
    ispy_b_static.prefix = "/bound/"
    ispy_b_static.input_string = file_path
    return ispy_b_static.get_response()


def map_download(request, file_path):
    """
    Download an event map by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = SiteObservation
    ispy_b_static.request = request
    ispy_b_static.permission_string = "target_id__project_id"

    substrings = file_path.split('.')[0].split('_')
    substring = [x for x in substrings if x in ['2fofc', 'fofc', 'event']]
    if not substring:
        file_extension = None
    else:
        file_extension = substring[0]

    # TODO: remove/add map_info (was used for hotspots but not currently used)

    exts = {'sigmaa_info': '2fofc', 'diff_info': 'fofc', 'event_info': 'event'}

    field_name = None

    for key in exts.keys():
        if exts[key] == file_extension:
            field_name = key

    ispy_b_static.field_name = field_name
    ispy_b_static.content_type = "application/x-pilot"
    ispy_b_static.prefix = "/maps/"
    ispy_b_static.input_string = file_path
    return ispy_b_static.get_response()


def metadata_download(request, file_path):
    """
    Download a metadata file by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = Target
    ispy_b_static.request = request
    ispy_b_static.permission_string = "project_id"
    ispy_b_static.field_name = "metadata"
    ispy_b_static.content_type = "application/x-pilot"
    ispy_b_static.prefix = "/metadata/"
    ispy_b_static.input_string = file_path
    return ispy_b_static.get_response()


def archive_download(request, file_path):
    """
    Download a protein by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = Target
    ispy_b_static.request = request
    ispy_b_static.permission_string = "project_id"
    ispy_b_static.field_name = "zip_archive"
    ispy_b_static.content_type = "application/zip"
    ispy_b_static.prefix = "/targets/"
    ispy_b_static.input_string = file_path
    ispy_b_static.file_format = 'raw'
    return ispy_b_static.get_response()
