from django.contrib.auth.models import User
from rest_framework.authtoken.models import Token
from django.core.exceptions import ObjectDoesNotExist
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rest_framework import viewsets
from django.http import HttpResponse
from viewer.models import Project


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


def draw_mol(smiles, height=200, width=200):
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
    drawer = Draw.MolDraw2DSVG(height, width)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return _transparentsvg(drawer.GetDrawingText().replace("svg:", ""))


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

    def get_proposals_for_user(self):
        user = self.request.user
        return self.get_proposals_for_user_from_django(user)

    def get_filter_dict(self, proposal_list):
        return {self.filter_permissions + "__title__in": proposal_list}


def get_params(smiles, request):
    height = None
    if "height" in request.GET:
        height = int(request.GET["height"])
    width = None
    if "width" in request.GET:
        width = int(request.GET["width"])
    return HttpResponse(draw_mol(smiles, width=width, height=height))


def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"].rstrip(".svg")
        return get_params(smiles, request)
    else:
        return HttpResponse("Please insert SMILES")
