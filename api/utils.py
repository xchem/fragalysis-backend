from django.contrib.auth.models import User
from django.http import Http404
from rest_framework.authtoken.models import Token
from django.core.exceptions import ObjectDoesNotExist
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from django.http import HttpResponse
from viewer.models import Project
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
import os
from rest_framework import viewsets


def get_conn():
    credentials = {
        "user": os.environ["ISPYB_USER"],
        "pw": os.environ["ISPYB_PASSWORD"],
        "host": os.environ["ISPYB_HOST"],
        "port": os.environ["ISPYB_PORT"],
        "db": "ispyb",
        "conn_inactivity": 360,
    }
    conn = Connector(**credentials)
    return conn


class ISpyBSafeQuerySet(viewsets.ReadOnlyModelViewSet):

    def get_queryset(self):
        """
        Optionally restricts the returned purchases to a given propsals
        """
        # The list of proposals this user can have
        proposal_list = self.get_proposals_for_user()
        # Add in the ones everyone has access to
        proposal_list.extend(self.get_open_proposals())
        # Must have a directy foreign key (project_id) for it to work
        filter_dict = self.get_filter_dict(proposal_list)
        return self.queryset.filter(**filter_dict).distinct()

    def get_open_proposals(self):
        """
        Returns the list of proposals anybody can access
        :return:
        """
        return ["OPEN"]

    def get_proposals_for_user_from_django(self, user):
        # Get the list of proposals for the user
        return list(
            Project.objects.filter(user_id=user.pk).values_list("title", flat=True)
        )

    def get_proposals_for_user_from_ispyb(self, user):
        with get_conn() as conn:
            core = conn.core
            rs = core.retrieve_sessions_for_person_login(user.username)
        prop_ids = [str(x["proposalId"]) + "-" + str(x["sessionNumber"]) for x in rs]
        return prop_ids

    def get_proposals_for_user(self):
        user = self.request.user
        get_from_ispyb = os.environ.get("ISPYB_FLAG", False)
        if get_from_ispyb:
            if user.is_authenticated:
                return self.get_proposals_for_user_from_ispyb(user)
            else:
                return []
        else:
            return self.get_proposals_for_user_from_django(user)

    def get_filter_dict(self, proposal_list):
        return {self.filter_permissions + "__title__in": proposal_list}


class ISpyBSafeStaticFiles:

    def get_queryset(self):
        query = ISpyBSafeQuerySet()
        query.request = self.request
        query.filter_permissions = self.permission_string
        query.queryset = self.model.objects.filter()
        queryset = query.get_queryset()
        return queryset

    def get_response(self):
        try:
            queryset = self.get_queryset()
            filter_dict = {self.field_name + "__endswith": self.input_string}
            object = queryset.get(**filter_dict)
            file_name = os.path.basename(str(getattr(object, self.field_name)))
            response = HttpResponse()
            response["Content-Type"] = self.content_type
            response["X-Accel-Redirect"] = self.prefix + file_name
            response["Content-Disposition"] = "attachment;filename=" + file_name
            return response
        except Exception:
            raise Http404
        return response


def get_token(request):
    """
    Get the authentication token for a givne request.
    Should just return an un-authenticated user token if nothing.
    :param request:
    :return:
    """
    try:
        user = User.objects.get(username=request.user)
        token, created = Token.objects.get_or_create(user=user)
        return token.key
    except ObjectDoesNotExist:
        return ""


def _transparentsvg(svg):
    """
    Give an SVG a white background
    :param svg:
    :return:
    """
    # Make the white background transparent
    tree = ET.fromstring(svg)
    rect = tree.find("rect")
    rect.set("style", rect.get("style").replace("#FFFFFF", "none"))
    # Recover some missing attributes for correct browser rendering
    tree.set("version", "1.1")
    tree.set("xmlns", "http://www.w3.org/2000/svg")
    tree.set("xmlns:rdkit", "http://www.rdkit.org/xml")
    tree.set("xmlns:xlink", "http://www.w3.org/1999/xlink")
    return '<?xml version="1.0" encoding="UTF-8"?>' + ET.tostring(tree).strip()


def draw_mol(smiles, height=200, width=200, img_type=None):
    """
    Draw a molecule from a smiles
    :param smiles: the SMILES to render
    :param height: the height in px
    :param width: the width in px
    :return: an SVG as a string of the inage
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "None Mol"
    AllChem.Compute2DCoords(mol)
    Chem.Kekulize(mol)
    if not height:
        height = 200
    if not width:
        width = 200
    if img_type == "png":
        options = DrawingOptions()
        # options.defaultColor = (1, 1, 1)
        options.bondLineWidth = 6
        img = Draw.MolToImage(mol, options=options)
        img = img.convert("RGBA")
        datas = img.getdata()
        newData = []
        for item in datas:
            if item[0] == 255 and item[1] == 255 and item[2] == 255:
                newData.append((255, 255, 255, 0))
            else:
                newData.append(item)
        img.putdata(newData)
        response = HttpResponse(content_type="image/png")
        img.save(response, "PNG")
        return response
    else:
        drawer = Draw.MolDraw2DSVG(height, width)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return _transparentsvg(drawer.GetDrawingText().replace("svg:", ""))


def get_params(smiles, request):
    height = None
    if "height" in request.GET:
        height = int(request.GET["height"])
    width = None
    if "width" in request.GET:
        width = int(request.GET["width"])
    img_type = request.GET.get("img_type", None)
    get_mol = draw_mol(smiles, width=width, height=height, img_type=img_type)
    if type(get_mol) == HttpResponse:
        return get_mol
    return HttpResponse(get_mol)


def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"].rstrip(".svg")
        return get_params(smiles, request)
    else:
        return HttpResponse("Please insert SMILES")
