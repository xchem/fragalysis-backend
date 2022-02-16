from celery import shared_task
from django.core.files.storage import default_storage
import pandas as pd
from rdkit.Chem import AllChem

from .validate import ValidateFile
from .createmodels import (
    createProjectModel,
    createBatchModel,
    createTargetModel,
    createMethodModel,
    createReactionModel,
    createProductModel,
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
        mculeids = []
        amounts = []
        project_id, project_name = createProjectModel(project_info)
        project_info["project_name"] = project_name

        grouped_targets = uploaded_df.groupby("batch-tag")
        
        for batch_tag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batch_tag=batch_tag, 
            )
            target_no = 1
            for target_smiles, target_mass in zip(
                group["targets"], group["amount-required-mg"]
            ):

                retrosynthesis_result = getManifoldretrosynthesis(target_smiles)
                routes = retrosynthesis_result["routes"]

                target_id = createTargetModel(
                    batch_id=batch_id,
                    smiles=target_smiles,
                    target_no=target_no,
                    target_mass=target_mass,
                )

                target_no += 1

                method_no = 1
                for route in routes:
                    no_steps = len(route["reactions"])

                    if no_steps > 0:
                        reactions = route["reactions"]

                        reactions_found = [
                            reaction for reaction in reactions if reaction["name"] in encoded_recipes
                        ]

                        if len(reactions_found) == no_steps:

                            method_id = createMethodModel(
                                target_id=target_id,
                                nosteps=no_steps,
                            )

                            product_no = 1
                            for reaction in reversed(reactions):
                                reaction_name = reaction["name"]
                                if reaction_name in encoded_recipes:
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
                                        target_no=target_no,
                                        method_no=method_no,
                                        product_no=product_no,
                                        product_smiles=product_smiles,
                                    )

                                    create_models = CreateEncodedActionModels(
                                        actions=actions,
                                        target_id=target_id,
                                        reaction_id=reaction_id,
                                        reactant_pair_smiles=reactant_smiles_ordered,
                                        reaction_name=reaction_name,
                                    )

                                    mculeids.append(create_models.mculeidlist)
                                    amounts.append(create_models.amountslist)

                                product_no += 1

                    method_no += 1

    # CreateMculeQuoteModel(mculeids=mculeids, amounts=amounts, project_id=project_id)

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
        mculeids = []
        amounts = []
        project_id, project_name = createProjectModel(project_info)
        project_info["project_name"] = project_name

        grouped_targets = uploaded_df.groupby("batch-tag")
        
        for batch_tag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batch_tag=batch_tag, 
            )

            target_no = 1
            product_no = 1
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
                    target_no=target_no,
                    target_mass=target_mass,
                )

                target_no += 1

                recipes = encoded_recipes[reaction_name]["recipes"]

                method_no = 1
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
                    target_no=target_no,
                    method_no=method_no,
                    product_no=product_no,
                    product_smiles=target_smiles,
                )

                create_models = CreateEncodedActionModels(
                    actions=actions,
                    target_id=target_id,
                    reaction_id=reaction_id,
                    reactant_pair_smiles=reactant_pair_smiles,
                    reaction_name=reaction_name,
                )

                mculeids.append(create_models.mculeidlist)
                amounts.append(create_models.amountslist)
                method_no += 1

    # CreateMculeQuoteModel(mculeids=mculeids, amounts=amounts, project_id=project_id)

    default_storage.delete(csv_fp)

    return validate_dict, validated, project_info