from api.utils import get_response
from viewer.models import Protein


def prot_download(request, file_path):
    model = Protein
    permission = "target_id__project_id"
    field_name = "pdb_info"
    content_type = "application/x-pilot"
    prefix = "/pdbs/"
    return get_response(model, permission, field_name, content_type, prefix, file_path)
