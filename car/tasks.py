from __future__ import annotations
from celery import shared_task, current_task
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from zipfile import ZipFile
import pandas as pd
from rdkit.Chem import AllChem

from car.models import (
    Batch,
    Target,
    Method,
    Reaction,
    CompoundOrder,
    OTScript,
    OTSession,
    OTBatchProtocol,   
)

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

from .opentrons.cartoot import CreateOTSession
from .opentrons.otwrite import otWrite


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

                    target_id = createTargetModel(
                        batch_id=batch_id,
                        smiles=target_smiles,
                        target_mass=target_mass,
                    )

                    retrosynthesis_result = getManifoldretrosynthesis(target_smiles)
                    routes = retrosynthesis_result["routes"]    
                    
                    if not routes:
                        continue
                    else:
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
                                        # Get lowest cost/leadtime/preferred vendor - TO DO!!!!
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


@shared_task(bind=True)
def createOTScript(batchids: list):
    """"
    Create otscripts and starting plates for a list of batch ids
    """ 
    protocol_summary = {}
    print(batchids)
    for batchid in batchids:
        allreactionquerysets = getBatchReactions(batchid=batchid)
        if allreactionquerysets:
            otbatchprotocol = OTBatchProtocol()
            otbatchprotocol.batch_id = Batch.objects.get(id=batchid)   
            otbatchprotocol.celery_task_id = current_task.request.id
            otbatchprotocol.save()

            maxsteps = findmaxlist(allreactionquerysets=allreactionquerysets)
            groupedreactionquerysets = groupReactions(
            allreactionquerysets=allreactionquerysets, maxsteps=maxsteps
            )
            for index, reactiongroup in enumerate(groupedreactionquerysets):
                if index == 0:
                    otsession = CreateOTSession(
                        otbatchprotocol=otbatchprotocol,
                        reactiongroupqueryset=reactiongroup,
                    )

                    otsessionobj = otsession.otsessionobj
                    alladdactionsquerysetflat = otsession.alladdactionquerysetflat
                    startingreactionplatequeryset = otsession.startingreactionplatequeryset

                    otWrite(
                        otsessionobj=otsessionobj,
                        alladdactionsquerysetflat=alladdactionsquerysetflat,
                    )
                if index > 0:
                    otsession = CreateOTSession(
                        otbatchprotocol=otbatchprotocol,
                        reactiongroupqueryset=reactiongroup,
                        inputplatequeryset=startingreactionplatequeryset,
                    )

                    otsessionobj = otsession.otsessionobj
                    alladdactionsquerysetflat = otsession.alladdactionquerysetflat

                    otWrite(
                        otsessionobj=otsessionobj,
                        alladdactionsquerysetflat=alladdactionsquerysetflat,
                    )    

            protocol_summary[batchid] = True
            
            batch_tag = getBatchTag(batchid=batchid)
            zipOTBatchProtocol(otbatchprotocol=otbatchprotocol, batch_tag=batch_tag)
            
        else:
            protocol_summary[batchid] = False

    return protocol_summary

def zipOTBatchProtocol(otbatchprotocol: Django_object, batch_tag: str):
    otsession_queryset = OTSession.objects.filter(otbatchprotocol_id=otbatchprotocol)
    for otsession_obj in otsession_queryset:
        compoundorder_queryset = CompoundOrder.object.filter(otsession_id=otsession_obj)
        otscript_queryset = OTScript.objects.filter(otsession_id=otsession_obj)

        compound_order_csvs = [ContentFile(compound_order_obj.ordercsv.read(), name=compound_order_obj.ordercsv.name) for compound_order_obj in compoundorder_queryset]
        otscripts = [ContentFile(otscript_obj.otscript.read(), name=otscript_obj.otscript.name) for otscript_obj in otscript_queryset]

        zipfile = ZipFile("otsession-{}-for-batch-{}.zip".format(otsession_obj.id, batch_tag), "w")
        [zipfile.write(compound_order_csv) for compound_order_csv in compound_order_csvs]
        [zipfile.write(otscript) for otscript in otscripts]
        zipfile_fn = default_storage.save(
        "otbatchprotocols/",
        ContentFile(zipfile),
    )
        otbatchprotocol.zipfile = zipfile_fn
        otbatchprotocol.save()
        zipfile.close()

def getTargets(batchid):
    targetqueryset = Target.objects.filter(batch_id=batchid).order_by("id")
    return targetqueryset


def getMethods(targetid):
    methodqueryset = Method.objects.filter(target_id=targetid).filter(otchem=True).order_by("id")
    return methodqueryset


def getReactions(methodid):
    reactionqueryset = Reaction.objects.filter(method_id=methodid).order_by("id")
    return reactionqueryset

def getBatchTag(batchid):
    batch_obj = Batch.objects.get(id=batchid)
    batch_tag = batch_obj.batch_tag
    return batch_tag

def getBatchReactions(batchid):
    targetqueryset = getTargets(batchid=batchid)
    if targetqueryset:
        allreactionquerysets = []
        for target in targetqueryset:
            methodqueryset = getMethods(targetid=target)
            if methodqueryset:
                for method in methodqueryset:
                    reactionqueryset = getReactions(methodid=method.id)
                    allreactionquerysets.append(reactionqueryset)
        return allreactionquerysets


def findnoallreactionsteps(allreactionquerysets: list):
    """ "
    Function to get all possible number of reactions steps for
    multiple methods
    """
    allnumberofsteps = set([len(x) for x in allreactionquerysets])
    return allnumberofsteps


def findmaxlist(allreactionquerysets: list):
    maxlength = max([len(i) for i in allreactionquerysets])
    return maxlength


def groupReactions(allreactionquerysets: list, maxsteps: int):
    """
    Groups reactionqueries into first reactions, second reactions and so on
    """
    groupedreactionquerysets = []
    for i in range(maxsteps):
        reactiongroup = [
            reactionqueryset[i]
            for reactionqueryset in allreactionquerysets
            if i <= len(reactionqueryset) - 1
        ]
        groupedreactionquerysets.append(reactiongroup)
    return groupedreactionquerysets
