import re
import xml.etree.ElementTree as ET
from typing import Optional, Tuple

from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponse
from frag.utils.network_utils import canon_input
from rdkit import Chem
from rdkit.Chem import AllChem, Atom, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rest_framework.authtoken.models import Token

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

# Target Access String (TAS) regular expression.
TAS_REGEX_RE: re.Pattern = re.compile(settings.TAS_REGEX)


def get_token(request):
    """
    Get the authentication token for a given request.
    Should just return an un-authenticated user token if nothing.
    :param request:
    :return:
    """
    try:
        user = User.objects.get(username=request.user)
        token, _ = Token.objects.get_or_create(user=user)
        return token.key
    except ObjectDoesNotExist:
        return ""


def validate_tas(tas: str) -> Tuple[bool, Optional[str]]:
    """Validate a Target Access String against the defined regular expression.
    If the match fails the TAS_ERROR_MSG is returned.
    """
    if TAS_REGEX_RE.match(tas):
        return True, None
    return False, settings.TAS_REGEX_ERROR_MSG


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


def highlight_diff(prb_mol, ref_mol, width, height):
    """
    Draw a molecule (prb_mol) with the differences from a reference model highlighted
    :param prb_mol: smiles of the probe molecule
    :param ref_mol: smiles of the reference molecule
    :param width: output image width
    :param height: output image height
    :return: svg string of the image
    """
    if not width:
        width = 200
    if not height:
        height = 200

    mols = [Chem.MolFromSmiles(prb_mol), Chem.MolFromSmiles(ref_mol)]
    for m in mols:
        Chem.Kekulize(m)
    match = Chem.rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, completeRingsOnly=True)
    match_mol = Chem.MolFromSmarts(match.smartsString)
    rdDepictor.Compute2DCoords(mols[0])
    unconserved = [
        i
        for i in range(mols[0].GetNumAtoms())
        if i not in mols[0].GetSubstructMatch(match_mol)
    ]

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mols[0], highlightAtoms=unconserved)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    return svg


def calc_bounds(conformer):
    """
    Calculate the X,Y bounds of a molecule, the minimum and maximum X and Y coordinates of any atom in the molecule.
    It is assumed that the molecule already has 2D coordinates.

    :param mol:
    :return: X and Y bounds as lists of length 2
    """
    x = [0, 0]
    y = [0, 0]
    num = conformer.GetOwningMol().GetNumAtoms()
    for i in range(num):
        pos = conformer.GetAtomPosition(i)
        if pos.x < x[0]:
            x[0] = pos.x
        if pos.x > x[1]:
            x[1] = pos.x
        if pos.y < y[0]:
            y[0] = pos.y
        if pos.y > y[1]:
            y[1] = pos.y
    return x, y


def calc_dims(x, y):
    """
    Given the X and Y bounds generated with calc_bounds() calculate the X and Y dimensions
    :param x:
    :param y:
    :return: The X and Y dimensions
    """
    xd = x[1] - x[0]
    yd = y[1] - y[0]
    return xd, yd


def draw_mol(
    smiles,
    height=49,
    width=150,
    bondWidth=1,
    scaling=1.0,
    img_type=None,
    highlightAtoms=None,
    atomcolors=None,
    highlightBonds=None,
    bondcolors=None,
    mol=None,
):
    """
    Generate a SVG image of a molecule specified as SMILES
    :param smiles: The molecuels as SMILES
    :param the height in px
    :param width: the width in px
    :param bondWidth:
    :param scaling:
    :param img_type: Ignored. Kept for now for backwards compatibility
    :param highlightAtoms:
    :param atomcolors:
    :param highlightBonds:
    :param bondcolors:
    :param mol: an optional mol to use instead of the smiles. Kept for now for backwards compatibility
    :return:
    """
    del img_type

    if highlightAtoms is None:
        highlightAtoms = []
    if atomcolors is None:
        atomcolors = []
    if highlightBonds is None:
        highlightBonds = []
    if bondcolors is None:
        bondcolors = {}

    if not mol:
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "None Mol"

    Chem.Kekulize(mol)
    if mol.GetNumConformers() == 0:
        idx = AllChem.Compute2DCoords(mol)
        conformer = mol.GetConformer(idx)
    else:
        conformer = mol.GetConformer(-1)  # the default conformer

    # Do some maths to determine the optimal font.
    # If you want to influence this use the scaling parameter.
    x, y = calc_bounds(conformer)
    dim_x, dim_y = calc_dims(x, y)
    # Protect scaling from Div0.
    # Scale factors are 0 if the corresponding dimension is not +ve, non=zero.
    scale_x = width / dim_x if dim_x > 0 else 0
    scale_y = height / dim_y if dim_y > 0 else 0
    scale = min(scale_x, scale_y)
    font = max(round(scale * scaling), 6)

    # Now we can generate the drawing
    rdMolDraw2D.PrepareMolForDrawing(mol, wedgeBonds=False)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    drawopt = drawer.drawOptions()
    drawopt.clearBackground = False
    drawopt.padding = 0.05
    drawopt.bondLineWidth = bondWidth
    drawopt.minFontSize = font
    drawopt.maxFontSize = font
    drawopt.additionalAtomLabelPadding = 0.15

    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlightAtoms,
        highlightAtomColors=atomcolors,
        highlightBonds=highlightBonds,
        highlightBondColors=bondcolors,
    )

    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def parse_vectors(vector_list):
    return [int(x) for x in vector_list.split(",")]


def parse_bool(input_string):
    if input_string.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif input_string.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise ValueError("Value not parsable")


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
    for i, _ in enumerate(spl_list):
        list_len = 4
        if i % list_len in [0, 1]:
            atom_ids.append(int(spl_list[i]))
        if i % list_len == 2:
            iso = int(spl_list[i])
        if i % list_len == 3:
            add_hs = parse_bool(spl_list[i])
            atom_id_1 = atom_ids[0]
            atom_id_2 = atom_ids[1]
            if add_hs:
                mol = AllChem.AddHs(mol)
                # Replace the H with the atom id in atom_ids[0], atom_ids[1] with *
                h_atoms = [x for x in mol.GetAtoms() if x.GetAtomicNum() == 1]
                atom_remove = [
                    x.GetIdx() for x in h_atoms if x.GetIdx() in [atom_id_1, atom_id_2]
                ][0]
                ed_mol = AllChem.EditableMol(mol)
                # Remove the other Hs
                ed_mol.ReplaceAtom(atom_remove, Atom(0))
                # Get a new editable molecule
                mol = ed_mol.GetMol()
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
                # Record the new Atom Ids
                atom_ids = [
                    [x.GetBonds()[0].GetBeginAtomIdx(), x.GetBonds()[0].GetEndAtomIdx()]
                    for x in mol.GetAtoms()
                    if x.GetAtomicNum() == 0
                ][0]
                atom_id_1 = atom_ids[0]
                atom_id_2 = atom_ids[1]
            bond = mol.GetBondBetweenAtoms(atom_id_1, atom_id_2)
            bond_ids.append(bond.GetIdx())
            bond_colours[bond.GetIdx()] = ISO_COLOUR_MAP[iso]
            atom_ids = []
    return bond_ids, bond_colours, mol


def parse_xenons(input_smi):
    mol = Chem.MolFromSmiles(input_smi)
    e_mol = AllChem.EditableMol(mol)
    xenons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 54]
    bond_ids = []
    bond_colours = {}
    for xe in xenons:
        bond_id = xe.GetBonds()[0].GetIdx()
        bond_ids.append(bond_id)
        if len(xenons) > 1:
            bond_colours[bond_id] = ISO_COLOUR_MAP[xe.GetIsotope()]
        else:
            bond_colours[bond_id] = ISO_COLOUR_MAP[101]
        e_mol.ReplaceAtom(xe.GetIdx(), Atom(0))
    return bond_ids, bond_colours, e_mol.GetMol()


def get_img_from_smiles(smiles, request):
    # try:
    smiles = canon_input(smiles)
    # except:
    #     smiles = ""
    mol = None
    bond_id_list = []
    highlightBondColors = {}
    height = int(request.GET.get("height", "128"))
    width = int(request.GET.get("width", "128"))
    if "atom_indices" in request.GET:
        mol = Chem.MolFromSmiles(smiles)
        bond_id_list, highlightBondColors, mol = parse_atom_ids(
            request.GET["atom_indices"], mol
        )
    if "Xe" in smiles:
        bond_id_list, highlightBondColors, mol = parse_xenons(smiles)
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


def get_highlighted_diffs(request):
    prb_smiles = request.GET['prb_smiles']
    ref_smiles = request.GET['ref_smiles']
    height = None
    width = None
    if "height" in request.GET:
        height = int(request.GET["height"])
    if "width" in request.GET:
        width = int(request.GET["width"])
    return HttpResponse(
        highlight_diff(
            prb_mol=prb_smiles, ref_mol=ref_smiles, height=height, width=width
        )
    )


def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"].rstrip(".svg")
        return get_img_from_smiles(smiles, request)
    else:
        return HttpResponse("Please insert SMILES")


def pretty_request(request, *, tag='', print_body=False):
    """A simple function to return a Django request as nicely formatted string."""
    headers = ''
    if request.headers:
        for key, value in request.headers.items():
            headers += f'{key}: {value}\n'

    tag_text = f'{tag}\n' if tag else ''

    user_text = 'User: '
    if request.user:
        user_text += str(request.user)
    else:
        user_text += '(-)'

    # Load the body but cater for problems, like
    #   django.http.request.RawPostDataException:
    #   You cannot access body after reading from request's data stream
    body = None
    if print_body:
        try:
            body = request.body
        except Exception:
            pass

    return (
        f'{tag_text}'
        '+ REQUEST BEGIN\n'
        f'{user_text}\n'
        f'{request.method} HTTP/1.1\n'
        f'{headers}\n\n'
        f'{body}\n'
        '- REQUEST END'
    )


def deployment_mode_is_production():
    return settings.DEPLOYMENT_MODE == "PRODUCTION"
