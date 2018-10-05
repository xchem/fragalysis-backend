import xml.etree.ElementTree as ET

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponse
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.Draw import rdMolDraw2D
from rest_framework.authtoken.models import Token
from frag.utils.network_utils import get_fragments, canon_input

ISO_COLOUR_MAP = {
    100: (1, 0, 0),
    101: (0, 1, 0),
    102: (0, 0, 1),
    103: (1, 0, 1),
    104: (1, 1, 0),
    105: (0, 1, 1),
    106: (0.5, 0.5, 0.5),
    107: (1, 0.5, 1),
}


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


def draw_mol(
    smiles,
    height=200,
    width=200,
    img_type=None,
    highlightAtoms=[],
    atomcolors=[],
    highlightBonds=[],
    bondcolors={},
    mol=None,
):
    """
    Draw a molecule from a smiles
    :param smiles: the SMILES to render
    :param height: the height in px
    :param width: the width in px
    :return: an SVG as a string of the inage
    """
    if mol is None:
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
        img = Draw.MolToImage(
            mol,
            options=options,
            highlightBonds=highlightBonds,
            highlightBondColors=bondcolors,
        )
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
        drawer = rdMolDraw2D.MolDraw2DSVG(height, width)
        drawopt = drawer.drawOptions()
        drawopt.clearBackground = False
        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlightAtoms,
            highlightAtomColors=atomcolors,
            highlightBonds=highlightBonds,
            highlightBondColors=bondcolors,
        )
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().replace("svg:", "")


def parse_vectors(vector_list):
    return [int(x) for x in vector_list.split(",")]


def parse_atom_ids(input_list, mol):
    """
    List of the form id,id,isotope,addHs
    e.g. 1,2,104,True
    :param input_list:
    :param mol:
    :return:
    """
    spl_list = input_list.split(",")
    bond_ids = []
    atom_ids = []
    bond_colours = {}
    for i, data in enumerate(spl_list):
        if i % 4 in [0, 1]:
            atom_ids.append(int(spl_list[i]))
        if i % 3 == 2:
            iso = int(spl_list[i])
        if i % 4 == 3:
            add_hs = bool(spl_list[i])
            if add_hs:
                mol = AllChem.AddHs(mol)
            bond = mol.GetBondBetweenAtoms(atom_ids[0], atom_ids[1])
            bond_ids.append(bond.GetIdx())
            bond_colours[bond.GetIdx] = ISO_COLOUR_MAP[iso]
            atom_ids = []
    return bond_ids, bond_colours, mol


def get_params(smiles, request):
    try:
        smiles = canon_input(smiles)
    except:
        smiles = ""
    height = None
    mol = None
    bond_id_list = []
    highlightBondColors = {}
    if "height" in request.GET:
        height = int(request.GET["height"])
    width = None
    if "width" in request.GET:
        width = int(request.GET["width"])
    if "isotopes" in request.GET:
        # Get the vector
        isotopes = parse_vectors(request.GET["isotopes"])
        mol = AllChem.AddHs(Chem.MolFromSmiles(smiles))
        fragments, frag_map = get_fragments(mol, get_index_iso_map=True)
        for iso in isotopes:
            if iso in frag_map:
                bond_ids = frag_map[iso]
                bond_id_list.extend(bond_ids)
                for bond_id in bond_ids:
                    highlightBondColors[bond_id] = ISO_COLOUR_MAP[iso]
    if "atom_indices" in request.GET:
        mol = Chem.MolFromSmiles(smiles)
        bond_ids, bond_colours, render_mol = parse_atom_ids(
            request.GET["atom_indices"], mol
        )
        mol = render_mol
        bond_id_list = bond_ids
        highlightBondColors = bond_colours
    img_type = request.GET.get("img_type", None)
    get_mol = draw_mol(
        smiles,
        width=width,
        height=height,
        img_type=img_type,
        highlightBonds=bond_id_list,
        mol=mol,
        bondcolors=highlightBondColors,
    )
    if type(get_mol) == HttpResponse:
        return get_mol
    return HttpResponse(get_mol)


def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"].rstrip(".svg")
        return get_params(smiles, request)
    else:
        return HttpResponse("Please insert SMILES")
