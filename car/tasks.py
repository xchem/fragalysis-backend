from celery import shared_task
from django.core.files.storage import default_storage
import pandas as pd
from rdkit.Chem import AllChem

from .validate import ValidateFile
from .createmodels import (
    createProjectModel,
    createBatchModel,
    createTargetModel,
    createCatalogEntryModel,
    createMethodModel,
    createReactionModel,
    createProductModel,
    createReactantModel,
    CreateEncodedActionModels,
    CreateMculeQuoteModel,
)

from .manifold.apicalls import getManifoldretrosynthesis
from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import getAddtionOrder


def delete_tmp_file(filepath):
    default_storage.delete(filepath)


# Celery validate task
@shared_task
def validateFileUpload(csv_fp, validate_type=None, project_info=None, validate_only=True):
    """Celery task to process validate the uploaded files for retrosynthesis planning.

    Parameters
    ----------
    csv_fp: str
        filepath of the uploaded csv file, which is saved to temporary storage by `viewer.views.UploadCSV`
    validate_type: validate different types of upload files
    project_info: dict
        dictionary of project details (name, email and project_name) that will be used to create the project model
        if the csv file is validated
    validate_only: boolean
        set to True to delete tmp file saved for validation and set to False to save tmp file and keep
        for uploading and creating model objects

    Returns
    -------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded csv file
    """

    validation = ValidateFile(csv_to_validate=csv_fp, validate_type=validate_type)

    if validate_only:
        default_storage.delete(csv_fp)
        csv_fp = None

    uploaded_dict = validation.df.to_dict("list")
    validated = validation.validated
    validate_dict = validation.validate_dict

    return (validate_dict, validated, csv_fp, project_info, uploaded_dict)


@shared_task
def uploadManifoldReaction(validate_output):

        validate_dict, validated, csv_fp, project_info, uploaded_dict = validate_output
        uploaded_df = pd.DataFrame(uploaded_dict)
   
        if not validated:
            default_storage.delete(csv_fp)
            return (validate_dict, validated, project_info)

        if validated:
            project_id, project_name = createProjectModel(project_info)
            project_info["project_name"] = project_name

            grouped_targets = uploaded_df.groupby("batch-tag")
            

            for batch_tag, group in grouped_targets:
                batch_id = createBatchModel(
                    project_id=project_id,
                    batch_tag=batch_tag, 
                )
  
                for target_smiles, target_mass in zip(
                    group["targets"], group["amount-required-mg"]
                ):

                    retrosynthesis_result = getManifoldretrosynthesis(target_smiles)
                    routes = retrosynthesis_result["routes"]    

                    target_id = createTargetModel(
                        batch_id=batch_id,
                        smiles=target_smiles,
                        target_mass=target_mass,
                    )
                    
                    first_route = routes[0]

                    # Check if target is in a catalogue and create catalog entries if it is
                    if first_route["molecules"][0]["isBuildingBlock"] == True:
                        catalog_entries = first_route["molecules"][0]["catalogEntries"]
                        for catalog_entry in catalog_entries:
                            createCatalogEntryModel(catalog_entry=catalog_entry, target_id=target_id)        

                    for route in routes[1:]:
                        no_steps = len(route["reactions"])
                        reactions = route["reactions"]
                        encoded_reactions_found = [
                            reaction for reaction in reactions if reaction["name"] in encoded_recipes
                        ]

                        if len(encoded_reactions_found) != no_steps:
                            method_id = createMethodModel(
                                target_id=target_id,
                                nosteps=no_steps,
                                otchem=False
                            )

                            for reaction in reversed(reactions):
                                reaction_name = reaction["name"]
                                reactant_smiles = reaction["reactantSmiles"]
                                product_smiles = reaction["productSmiles"]
                                reaction_smarts = AllChem.ReactionFromSmarts(
                                    "{}>>{}".format(
                                        ".".join(reactant_smiles), product_smiles
                                    ),
                                    useSmiles=True,
                                )

                                reaction_id = createReactionModel(
                                    method_id=method_id,
                                    reaction_class=reaction_name,
                                    reaction_smarts=reaction_smarts,
                                )

                                createProductModel(
                                    reaction_id=reaction_id,
                                    project_name=project_name,
                                    batch_tag=batch_tag,
                                    product_smiles=product_smiles,
                                )

                                for reactant_smi in reactant_smiles:
                                    reactant_id = createReactantModel(
                                                    reaction_id=reaction_id,
                                                    reactant_smiles=reactant_smi,
                                                    )
                                    catalog_entries = [molecule["catalogEntries"] for molecule in route["molecules"] if molecule["smiles"] == reactant_smi][0]
                                    for catalog_entry in catalog_entries:
                                        createCatalogEntryModel(
                                            catalog_entry=catalog_entry, 
                                            reactant_id=reactant_id
                                        )

                        if len(encoded_reactions_found) == no_steps:
                            method_id = createMethodModel(
                                target_id=target_id,
                                nosteps=no_steps,
                                otchem=True
                            )

                            for reaction in reversed(reactions):
                                reaction_name = reaction["name"]
                                recipes = encoded_recipes[reaction_name]["recipes"]
                                recipe_rxn_smarts = encoded_recipes[reaction_name]["reactionSMARTS"]
                                reactant_smiles = reaction["reactantSmiles"]
                                product_smiles = reaction["productSmiles"]

                                if len(reactant_smiles) == 1:
                                    actions = recipes["Intramolecular"]["actions"]
                                    stir_action = [
                                        action for action in actions if action["name"] == "stir"
                                    ][0]
                                    reaction_temperature = stir_action["content"]["temperature"][
                                        "value"
                                    ]
                                    reactant_smiles_ordered = reactant_smiles
                                else:
                                    actions = recipes["Standard"]["actions"]
                                    stir_action = [
                                        action for action in actions if action["name"] == "stir"
                                    ][0]
                                    reaction_temperature = stir_action["content"]["temperature"][
                                        "value"
                                    ]
                                    reactant_smiles_ordered = getAddtionOrder(
                                        product_smi=product_smiles,
                                        reactant_SMILES=reactant_smiles,
                                        reaction_SMARTS=recipe_rxn_smarts,
                                    )
                                    if not reactant_smiles_ordered:
                                        continue

                                reaction_smarts = AllChem.ReactionFromSmarts(
                                    "{}>>{}".format(
                                        ".".join(reactant_smiles_ordered), product_smiles
                                    ),
                                    useSmiles=True,
                                )

                                reaction_id = createReactionModel(
                                    method_id=method_id,
                                    reaction_class=reaction_name,
                                    reaction_temperature=reaction_temperature,
                                    reaction_smarts=reaction_smarts,
                                )

                                createProductModel(
                                    reaction_id=reaction_id,
                                    project_name=project_name,
                                    batch_tag=batch_tag,
                                    product_smiles=product_smiles,
                                )

                                for reactant_smi in reactant_smiles_ordered:
                                    reactant_id = createReactantModel(
                                                    reaction_id=reaction_id,
                                                    reactant_smiles=reactant_smi,
                                                    )
                                    
                                    catalog_entries = [molecule["catalogEntries"] for molecule in route["molecules"] if molecule["smiles"] == reactant_smi][0]
                                    for catalog_entry in catalog_entries:
                                        createCatalogEntryModel(
                                            catalog_entry=catalog_entry, 
                                            reactant_id=reactant_id
                                        )


                                CreateEncodedActionModels(
                                    actions=actions,
                                    target_id=target_id,
                                    reaction_id=reaction_id,
                                    reactant_pair_smiles=reactant_smiles_ordered,
                                    reaction_name=reaction_name,
                                )

            default_storage.delete(csv_fp)

            return validate_dict, validated, project_info


@shared_task
def uploadCustomReaction(validate_output):

    validate_dict, validated, csv_fp, project_info, uploaded_dict = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)

    if not validated:
        default_storage.delete(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        project_id, project_name = createProjectModel(project_info)
        project_info["project_name"] = project_name

        grouped_targets = uploaded_df.groupby("batch-tag")
        
        for batch_tag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batch_tag=batch_tag, 
            )

            for reactant_pair_smiles, reaction_name, target_smiles, target_mass in zip(
                group["reactant-pair-smiles"],
                group["reaction-name"],
                group["target-smiles"],
                group["amount-required-mg"],
            ):
                reaction_smarts = AllChem.ReactionFromSmarts(
                    "{}>>{}".format(".".join(reactant_pair_smiles), target_smiles),
                    useSmiles=True,
                )

                target_id = createTargetModel(
                    batch_id=batch_id,
                    smiles=target_smiles,
                    target_mass=target_mass,
                )

                recipes = encoded_recipes[reaction_name]["recipes"]

                actions = recipes["Standard"]["actions"]
                stir_action = [action for action in actions if action["name"] == "stir"][0]
                reaction_temperature = stir_action["content"]["temperature"]["value"]

                method_id = createMethodModel(
                    target_id=target_id,
                    nosteps=1,
                )

                reaction_id = createReactionModel(
                    method_id=method_id,
                    reaction_class=reaction_name,
                    reaction_temperature=reaction_temperature,
                    reaction_smarts=reaction_smarts,
                )

                createProductModel(
                    reaction_id=reaction_id,
                    project_name=project_name,
                    batch_tag=batch_tag,
                    product_smiles=target_smiles,
                )

                CreateEncodedActionModels(
                    actions=actions,
                    target_id=target_id,
                    reaction_id=reaction_id,
                    reactant_pair_smiles=reactant_pair_smiles,
                    reaction_name=reaction_name,
                )

    default_storage.delete(csv_fp)

    return validate_dict, validated, project_info