from django.db.models import QuerySet
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
import itertools
import re
import inspect

import logging

from car.models import Batch, Method, Product, Reaction

from car.recipebuilder.encodedrecipes import encoded_recipes

logger = logging.getLogger(__name__)


# def getReactionInfo(batch_ids: list[int]) -> list:
#     """Get the reaction info for a batch as list of targets,reactant SMILES,
#     reaction class (Amidation, Suzuki etc)
#     """
#     reaction_info = {}
#     for batch_id in batch_ids:
#         batchobj = Batch.objects.get(id=batch_id)
#         targetqueryset = batchobj.targets.all()
#         target_smiles = [target_obj.smiles for target_obj in targetqueryset]
#         for targetobj in targetqueryset:
#             reaction_info["target_smiles"] = targetobj.smiles

#         methods = [target.methods.all() for target in targetqueryset]
#         method_sublist = [item for sublist in methods for item in sublist]
#         reactions = [method.reactions.all() for method in method_sublist]
#         reaction_sublist = [item for sublist in reactions for item in sublist]
#         reactants = [reaction.reactants.all() for reaction in reaction_sublist]
#         reactants_sublist = [item for sublist in reactants for item in sublist]
#         reactants_batch_to_buy = list(
#             set(
#                 [
#                     reactant.smiles
#                     for reactant in reactants_sublist
#                     if reactant.previousreactionproduct == False
#                 ]
#             )
#         )
#         reactants_to_buy = reactants_to_buy + reactants_batch_to_buy


# Need to convert into function to retrieve all reactants that need to be purchased plus API endpoint!
def getReactantsToBuy(batch_ids: list[int]) -> list:
    """Finds the reactnats that need to be bought to execute a batch/batches
    synthesis. Finds recatants that are not made in previous method's reactions

    Parameters
    ----------
    batch_ids: list[int]
        The batch ids to search for reactants to buy to complete the synthesis

    Returns
    -------
    reactants_to_buy: list
        The SMILES of the reactants that need to be bought. Excludes reactants
        made in previous reaction steps
    """
    reactants_to_buy = []
    for batch_id in batch_ids:
        batchobj = Batch.objects.get(id=batch_id)
        targetqueryset = batchobj.targets.all()
        methods = [target.methods.all() for target in targetqueryset]
        method_sublist = [item for sublist in methods for item in sublist]
        reactions = [method.reactions.all() for method in method_sublist]
        reaction_sublist = [item for sublist in reactions for item in sublist]
        reactants = [reaction.reactants.all() for reaction in reaction_sublist]
        reactants_sublist = [item for sublist in reactants for item in sublist]
        reactants_batch_to_buy = list(
            set(
                [
                    reactant.smiles
                    for reactant in reactants_sublist
                    if reactant.previousreactionproduct == False
                ]
            )
        )
        reactants_to_buy = reactants_to_buy + reactants_batch_to_buy
    return list(set(reactants_to_buy))


def getProduct(reaction_id: int) -> Product:
    """Get product object

    Parameters
    ----------
    reaction_id: int
        The reaction id to search for a matching product

    Returns
    -------
    productobj: Product
        The product Django model object
    """
    productobj = Product.objects.get(reaction_id=reaction_id)
    return productobj


def getReaction(reaction_id: int) -> Reaction:
    """Get reaction object

    Parameters
    ----------
    reaction_id: int
        The reaction id to search for a reaction

    Returns
    -------
    reactionobj: Reaction
        The reaction Django model object
    """
    reactionobj = Reaction.objects.get(id=reaction_id)
    return reactionobj


def getReactionQuerySet(
    reaction_ids: list = None, method_id: int = None
) -> QuerySet[Reaction]:
    """Get a  synthesis methods reactions

    Parameters
    ----------
    reaction_id: int or Reaction
        The reaction ids to find reactions for
    method_id: int
        The optional synthesis method's id to get reactions for

    Returns
    -------
    reactionqueryset: QuerySet[Reaction]
        The reactions of a synthesis method
    """
    if reaction_ids:
        reactionqueryset = Reaction.objects.filter(id__in=reaction_ids).order_by("id")
    if method_id:
        reactionqueryset = Reaction.objects.filter(method_id=method_id).order_by("id")
    return reactionqueryset


def checkProceedingReactions(reaction_id: int) -> QuerySet[Reaction]:
    """Checks if there are any reactions that proceed the reaction

    Parameters
    ----------
    reaction_id: int
        The reaction id of the Django model object to search for
        all relative proceeding reactions objects

    Returns
    -------
    proceedingreactionqueryset: QuerySet[Reaction]
        Returns the reactions that proceed the reaction
    """
    reactionobj = getReaction(reaction_id=reaction_id)
    proceedingreactionqueryset = Method.objects.get(
        id=reactionobj.method_id.id
    ).reactions.filter(id__gt=reaction_id)
    return proceedingreactionqueryset


def getReactionYields(reactionclasslist: list) -> list[int]:
    """Gets the reaction yields

    Parameters
    ----------
    reactionclasslist: list
        The reaction classes to find yields for

    Returns
    -------
    reactionyields: list[float]
        Returns the reaction yields eg. 0.80
    """
    reactionyields = [
        (encoded_recipes[reactionclass]["recipes"]["standard"]["yield"] / 100)
        for reactionclass in reactionclasslist
    ]
    return reactionyields


def checkPreviousReactionProducts(reaction_id: int, smiles: str) -> QuerySet[Reaction]:
    """Checks if any previous reactions had a product matching the smiles

    Parameters
    ----------
    reaction_id: int
        The reaction id of the Django model object to search for
        all relative previous reactions objects. The previosu reactions may
        have products that are this reaction's reactant input
    smiles: str
        The SMILES of the reaction's reactant and previous reaction products

    Returns
    -------
    previousreactionqueryset: QuerySet[Reaction]
        Returns the reactions that yiled products that match the SMILES searched
    """
    reactionobj = getReaction(reaction_id=reaction_id)
    previousreactionqueryset = Method.objects.get(
        id=reactionobj.method_id.id
    ).reactions.filter(id__lt=reaction_id, products__smiles=smiles)
    return previousreactionqueryset


def getMWs(smiles: list[str]) -> list[float]:
    """Gets the molecular weights of a list of compounds SMILES

    Parameters
    ----------
    smiles: list[str]
        The SMILES to calculate molecular weights for

    Returns
    -------
    MWs: list[float]
        The list of molecular weights
    """

    MWs = [Descriptors.ExactMolWt(Chem.MolFromSmiles(smi)) for smi in smiles]
    return MWs


def calculateProductMols(target_mass: float, target_SMILES: str) -> object:
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


def combiChem(reactant_1_SMILES: list, reactant_2_SMILES: list) -> list:
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
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
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
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
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
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
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
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(
            "Pubchempy could not convert SMILES to a IUPAC name with error {}".format(e)
        )
        return None
