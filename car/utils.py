from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
import itertools
import re


def calculateproductmols(target_mass: float, target_SMILES: str) -> object:
    """Function to calculate product mols of reaction using a target mass

    Parameters
    ----------
    target_mass: float
        The target mass (mg) of the product
    target_SMILES: str
        The SMILES of the product

    Returns
    -------
    product_moles: rdkit mol object
        The product mols
    """
    target_MW = Descriptors.ExactMolWt(Chem.MolFromSmiles(target_SMILES))
    target_mass = target_mass / 1e3
    product_mols = target_mass / target_MW
    return product_mols


def canonSmiles(smiles: str) -> str:
    """Function to canonicalise SMILES

    Parameters
    ----------
    smiles: str
        The SMILES to be canonicalised

    Returns
    -------
    canon_smiles: str
        The canonicalised SMILES
    status: bool
        Returns False if the input smiles could not be canonicalised
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        canon_smiles = Chem.MolToSmiles(mol)
        return canon_smiles
    else:
        return False


def combichem(reactant_1_SMILES: list, reactant_2_SMILES: list) -> list:
    """Gets all possible combinations between two uneven lists of
       reactants

    Parameters
    ----------
    reactant_1_SMILES: list
        The list of reactant one smiles
    reactant_2_SMILES: list
        The second list of reactant two smiles

    Returns
    -------
    all_possible_combinations: list
        All possible reactant combinations possible
        between reactat 1 and reactant two lists
        as a list of tuples
    """
    all_possible_combinations = list(
        itertools.product(reactant_1_SMILES, reactant_2_SMILES)
    )

    return all_possible_combinations


def createSVGString(smiles: str) -> str:
    """Function that creates a SVG image string from smiles string

    Parameters
    ----------
    smiles: string
        The SMILES to create an SVG image string from

    Returns
    -------
    svg_string: string
        The SVG image string
    """
    mol = Chem.MolFromSmiles(smiles)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(100, 50)
    drawer.SetFontSize(8)
    drawer.SetLineWidth(1)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg_string = drawer.GetDrawingText()
    return svg_string


def createReactionSVGString(smarts: str) -> str:
    """Function that creates a SVG image string from smarts string

    Parameters
    ----------
    smarts: string
        The SMARTS reaction pattern to create an SVG image
        string from

    Returns
    -------
    svg_string: string
        The SVG image string of the SMARTS pattern
    """
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(900, 200)
    drawer.DrawReaction(smarts)
    drawer.FinishDrawing()
    svg_string = drawer.GetDrawingText()
    return svg_string


def getAddtionOrder(
    product_smi: str, reactant_SMILES: tuple, reaction_SMARTS: str
) -> list:
    """Gets reactant pair addition order as SMILES that yields the expected
       prodcut via the reaction SMARTS pattern

    Parameters
    ----------
    product_smi: str
        The product SMILES
    reactant_SMILES_pair: tuple
        The reactant SMILES pair for a reaction
    reaction_SMARTS: str
        The reaction SMARTS pattern

    Returns
    -------
    reactant_SMILES_pair: list
        The list of ordered reactant smiles
    status: None
        None if no order can create the input product
    """
    rxn = AllChem.ReactionFromSmarts(reaction_SMARTS)
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactant_SMILES]

    for reactant_permutation in list(itertools.permutations(reactant_mols)):
        try:
            products = rxn.RunReactants(reactant_permutation)
            product_mols = [product[0] for product in products]
            if not product_mols:
                continue  # reactants were in wrong order so no product
        except Exception as e:
            print(e)
            print(reactant_permutation)
            continue
        product_smis = [Chem.MolToSmiles(m) for m in product_mols if m is not None]
        if product_smi in product_smis:
            ordered_smis = [Chem.MolToSmiles(m) for m in reactant_permutation]
    if "ordered_smis" in locals():
        return ordered_smis
    else:
        print(reaction_SMARTS)
        print(reactant_SMILES)
        return None


def checkReactantSMARTS(reactant_SMILES: tuple, reaction_SMARTS: str) -> list:
    """Checks if reactant pair can produce a product

    Parameters
    ----------
    reactant_SMILES_pair: tuple
        The pair of reactant smiles to check
    reaction_SMARTS: str
        The reaction SMARTS pattern used to check the reactant SMILES

    Returns
    -------
    product_mols: list
        The list of product mols formed between reactant SMILES from SMARTS pattern
    status: None
        Returns None if no product mols are formed
    """
    rxn = AllChem.ReactionFromSmarts(reaction_SMARTS)
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactant_SMILES]

    for reactant_permutation in list(itertools.permutations(reactant_mols)):
        try:
            products = rxn.RunReactants(reactant_permutation)
            product_mols = [product[0] for product in products]
            if product_mols:
                break
            if not product_mols:
                continue  # reactants were in wrong order so no product
        except Exception as e:
            print(e)
            print(reactant_permutation)
            continue
    if "product_mols" in locals():
        return product_mols
    else:
        print(reaction_SMARTS)
        print(reactant_SMILES)
        return None


def getPubChemCAS(compound: object) -> str:
    """Get CAS identifier for PubChem compound synonyms

    Parameters
    ----------
    compound: object
        A PuBChem compound object

    Returns
    -------
    cas: str
        The CAS id of the compound
    """
    synonyms = compound.synonyms
    if synonyms:
        for syn in synonyms:
            match = re.match("(\d{1,7}-\d{1,2}-\d)", syn)
            if match:
                cas = match.group(1)
                return cas


def getPubChemCompound(smiles: str) -> object:
    """Searches PubChem for compound using SMILES

    Parameters
    ----------
    smiles: str
        The SMILES of the compound to search the PubChem DB for

    Returns
    -------
    compound: object
        The PuBChem compound class object
    status: None
        Returns None if no compound is found or an error occurs
    """
    try:
        compound = pcp.get_compounds(smiles, "smiles")[0]
        if not compound.cid:
            return None
        else:
            return compound
    except Exception as e:
        print(
            "Pubchempy could not retrieve compound entry for input SMILES: {} with error {}".format(
                smiles, e
            )
        )
        return None


def getChemicalName(smiles: str) -> str:
    """Searches PubChem for compound using SMILES

    Parameters
    ----------
    smiles: str
        The SMILES of the compound to search the PubChem DB for
        it's IUPAC name

    Returns
    -------
    name: str
        The IUPAC name of the compound
    status: None
        Returns None if no compound IUPAC name is found or if an error
        occurs
    """
    try:
        name = pcp.get_compounds(smiles, "smiles")[0].iupac_name
        if not name:
            return None
        else:
            return name
    except Exception as e:
        print(
            "Pubchempy could not convert SMILES to a IUPAC name with error {}".format(e)
        )
        return None
