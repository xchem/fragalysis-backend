from django.contrib.auth.models import User
from rest_framework.authtoken.models import Token
from django.core.exceptions import ObjectDoesNotExist
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import AllChem,Draw

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
    rect = tree.find('rect')
    rect.set('style', rect.get('style').replace('#FFFFFF', 'none'))
    # Recover some missing attributes for correct browser rendering
    tree.set('version', '1.1')
    tree.set('xmlns', 'http://www.w3.org/2000/svg')
    tree.set('xmlns:rdkit', 'http://www.rdkit.org/xml')
    tree.set('xmlns:xlink', 'http://www.w3.org/1999/xlink')
    return '<?xml version="1.0" encoding="UTF-8"?>' + ET.tostring(tree).strip()

def draw_mol(smiles,height=200,width=200):
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
        height=200
    if not width:
        width=200
    drawer = Draw.MolDraw2DSVG(height,width)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return _transparentsvg(drawer.GetDrawingText().replace('svg:',''))