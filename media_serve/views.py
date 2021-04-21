from api.security import ISpyBSafeStaticFiles
from viewer.models import Protein, Target


def prot_download(request, file_path):
    """
    Download a protein by nginx redirect
    :param request: the initial request
    :param file_path: the file path we're getting from the static
    :return: the response (a redirect to nginx internal)
    """
    ispy_b_static = ISpyBSafeStaticFiles()
    ispy_b_static.model = Protein
    ispy_b_static.request = request
    ispy_b_static.permission_string = "target_id__project_id"
    ispy_b_static.field_name = "pdb_info"
    ispy_b_static.content_type = "application/x-pilot"
    ispy_b_static.prefix = "/pdbs/"
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
    ispy_b_static.model = Protein
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
    ispy_b_static.model = Protein
    ispy_b_static.request = request
    ispy_b_static.permission_string = "target_id__project_id"

    substrings = file_path.split('_')
    substring = [x for x in substrings if x in ['2fofc', 'fofc', 'event']]
    if not substring:
        file_extension = None
    else:
        file_extension = substring[0]

    # TODO: remove/add map_info (was used for hotspots but not currently used)

    exts = {'sigmaa_info': '2fofc',
            'diff_info': 'fofc',
            'event_info': 'event'}

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