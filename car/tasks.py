"""CAR's celery tasks"""
from __future__ import annotations
from celery import shared_task, current_task
from django.conf import settings
from django.db.models import QuerySet
from django.db.models import Q, Max

from zipfile import ZipFile
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os

from car.models import (
    ActionSession,
    Batch,
    SolventPrep,
    Target,
    Method,
    Reaction,
    CompoundOrder,
    OTProject,
    OTBatchProtocol,
    OTScript,
    OTSession,
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
)

from .manifold.apicalls import (
    getExactSearch,
    getManifoldRetrosynthesisBatch,
)
from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import checkPreviousReactionProducts, getAddtionOrder, canonSmiles

from .opentrons.otsession import CreateOTSession
from .opentrons.otwrite import OTWrite


def delete_tmp_file(filepath):
    os.remove(filepath)


@shared_task
def validateFileUpload(
    csv_fp, validate_type: str = None, project_info=None, validate_only=True
):
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
        delete_tmp_file(csv_fp)
        csv_fp = None

    uploaded_dict = validation.df.to_dict("list")
    validated = validation.validated
    validate_dict = validation.validate_dict

    return (
        validate_type,
        validate_dict,
        validated,
        project_info,
        csv_fp,
        uploaded_dict,
    )


@shared_task
def uploadManifoldReaction(validate_output):

    _, validate_dict, validated, project_info, csv_fp, uploaded_dict = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)

    if not validated:
        delete_tmp_file(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        project_id = createProjectModel(project_info)
        project_info["project_id"] = project_id

        grouped_targets = uploaded_df.groupby("batch-tag")

        for batchtag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batchtag=batchtag,
            )
            target_smiles = list(group["targets"])
            target_amounts = list(group["amount-required-mg"])
            for i in range(
                0, len(target_smiles), 10
            ):  # Manifold can do 10 smiles in one batch query
                smiles = target_smiles[i : i + 10]
                amounts = target_amounts[i : i + 10]
                retrosynthesis_results = getManifoldRetrosynthesisBatch(smiles=smiles)
                if "results" in retrosynthesis_results:
                    batchrouteresults = retrosynthesis_results["results"]

                    for smiles, mass, routeresult in zip(
                        smiles, amounts, batchrouteresults
                    ):
                        if "routes" in routeresult:
                            routes = routeresult["routes"]
                            if routes:
                                target_id = createTargetModel(
                                    batch_id=batch_id,
                                    smiles=smiles,
                                    mass=mass,
                                )
                                first_route = routes[0]

                                if first_route["molecules"][0]["isBuildingBlock"]:
                                    catalog_entries = first_route["molecules"][0][
                                        "catalogEntries"
                                    ]
                                    for catalog_entry in catalog_entries:
                                        createCatalogEntryModel(
                                            catalog_entry=catalog_entry,
                                            target_id=target_id,
                                        )

                                for route in routes[1:]:
                                    no_steps = len(route["reactions"])
                                    reactions = route["reactions"]
                                    encoded_reactions_found = [
                                        reaction
                                        for reaction in reactions
                                        if reaction["name"] in encoded_recipes
                                    ]

                                    if len(encoded_reactions_found) != no_steps:
                                        method_id = createMethodModel(
                                            target_id=target_id,
                                            nosteps=no_steps,
                                            otchem=False,
                                        )

                                        for index, reaction in enumerate(
                                            reversed(reactions)
                                        ):
                                            reaction_name = reaction["name"]
                                            reactant_smiles = reaction["reactantSmiles"]
                                            product_smiles = reaction["productSmiles"]
                                            reaction_smarts = (
                                                AllChem.ReactionFromSmarts(
                                                    "{}>>{}".format(
                                                        ".".join(reactant_smiles),
                                                        product_smiles,
                                                    ),
                                                    useSmiles=True,
                                                )
                                            )
                                            if len(reactant_smiles) == 1:
                                                intramolecular = True
                                            if len(reactant_smiles) == 2:
                                                intramolecular = False

                                            reaction_id = createReactionModel(
                                                method_id=method_id,
                                                reaction_class=reaction_name,
                                                reaction_number=index + 1,
                                                intramolecular=intramolecular,
                                                reaction_smarts=reaction_smarts,
                                            )

                                            createProductModel(
                                                reaction_id=reaction_id,
                                                product_smiles=product_smiles,
                                            )

                                            for reactant_smi in reactant_smiles:
                                                previousreactionqueryset = (
                                                    checkPreviousReactionProducts(
                                                        reaction_id=reaction_id,
                                                        smiles=reactant_smi,
                                                    )
                                                )
                                                if previousreactionqueryset:
                                                    previous_reaction_product = True
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    createCatalogEntryModel(
                                                        reactant_id=reactant_id,
                                                        previous_reaction_product=True,
                                                    )

                                                if not previousreactionqueryset:
                                                    previous_reaction_product = False
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    catalog_entries = [
                                                        molecule["catalogEntries"]
                                                        for molecule in route[
                                                            "molecules"
                                                        ]
                                                        if molecule["smiles"]
                                                        == reactant_smi
                                                    ][0]
                                                    for (
                                                        catalog_entry
                                                    ) in catalog_entries:
                                                        createCatalogEntryModel(
                                                            catalog_entry=catalog_entry,
                                                            reactant_id=reactant_id,
                                                        )

                                    if len(encoded_reactions_found) == no_steps:
                                        method_id = createMethodModel(
                                            target_id=target_id,
                                            nosteps=no_steps,
                                            otchem=True,
                                        )

                                        for index, reaction in enumerate(
                                            reversed(reactions)
                                        ):
                                            reaction_name = reaction["name"]
                                            reactant_smiles = reaction["reactantSmiles"]
                                            product_smiles = reaction["productSmiles"]
                                            reaction_temperature = [
                                                actionsession["actions"][0]["content"][
                                                    "temperature"
                                                ]["value"]
                                                for actionsession in encoded_recipes[
                                                    reaction_name
                                                ]["recipes"]["standard"][
                                                    "actionsessions"
                                                ]
                                                if actionsession["type"] == "stir"
                                            ][0]
                                            intramolecular_possible = encoded_recipes[
                                                reaction_name
                                            ]["intramolecular"]

                                            recipe_rxn_smarts = encoded_recipes[
                                                reaction_name
                                            ]["recipes"]["standard"]["reactionSMARTS"]

                                            if (
                                                len(reactant_smiles) == 1
                                                and intramolecular_possible
                                            ):
                                                reactant_smiles_ordered = (
                                                    reactant_smiles
                                                )
                                                intramolecular = True
                                            else:
                                                reactant_smiles_ordered = getAddtionOrder(
                                                    product_smi=product_smiles,
                                                    reactant_SMILES=reactant_smiles,
                                                    reaction_SMARTS=recipe_rxn_smarts,
                                                )
                                                intramolecular = False
                                                if not reactant_smiles_ordered:
                                                    continue

                                            reaction_smarts = (
                                                AllChem.ReactionFromSmarts(
                                                    "{}>>{}".format(
                                                        ".".join(
                                                            reactant_smiles_ordered
                                                        ),
                                                        product_smiles,
                                                    ),
                                                    useSmiles=True,
                                                )
                                            )

                                            reaction_id = createReactionModel(
                                                method_id=method_id,
                                                reaction_class=reaction_name,
                                                reaction_number=index + 1,
                                                intramolecular=intramolecular,
                                                recipe_type="standard",  # Need to change when we get multiple recipes runnning!
                                                reaction_temperature=reaction_temperature,
                                                reaction_smarts=reaction_smarts,
                                            )

                                            createProductModel(
                                                reaction_id=reaction_id,
                                                product_smiles=product_smiles,
                                            )

                                            for reactant_smi in reactant_smiles_ordered:
                                                previousreactionqueryset = (
                                                    checkPreviousReactionProducts(
                                                        reaction_id=reaction_id,
                                                        smiles=reactant_smi,
                                                    )
                                                )
                                                if previousreactionqueryset:
                                                    previous_reaction_product = True
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    createCatalogEntryModel(
                                                        reactant_id=reactant_id,
                                                        previous_reaction_product=True,
                                                    )

                                                if not previousreactionqueryset:
                                                    previous_reaction_product = False
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    catalog_entries = [
                                                        molecule["catalogEntries"]
                                                        for molecule in route[
                                                            "molecules"
                                                        ]
                                                        if molecule["smiles"]
                                                        == reactant_smi
                                                    ][0]
                                                    for (
                                                        catalog_entry
                                                    ) in catalog_entries:
                                                        createCatalogEntryModel(
                                                            catalog_entry=catalog_entry,
                                                            reactant_id=reactant_id,
                                                        )

        delete_tmp_file(csv_fp)

        return validate_dict, validated, project_info


@shared_task
def uploadCustomReaction(validate_output):

    (
        _,
        validate_dict,
        validated,
        project_info,
        csv_fp,
        uploaded_dict,
    ) = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)

    if not validated:
        delete_tmp_file(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        project_id = createProjectModel(project_info)
        project_info["project_id"] = project_id

        grouped_targets = uploaded_df.groupby("batch-tag")
        for batchtag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batchtag=batchtag,
            )

            for (
                target_smiles,
                amount_required,
                no_steps,
                reactant_pair_smiles_tuples,
                reaction_name_tuples,
                reaction_product_smiles_tuples,
            ) in zip(
                group["target-smiles"],
                group["amount-required-mg"],
                group["no-steps"],
                group["reactant-pair-smiles"],
                group["reaction-name"],
                group["product-smiles"],
            ):
                target_id = createTargetModel(
                    batch_id=batch_id,
                    smiles=target_smiles,
                    mass=amount_required,
                )
                method_id = createMethodModel(
                    target_id=target_id,
                    nosteps=no_steps,
                    otchem=True,
                )

                for reactant_pair_smiles, reaction_name, reaction_product_smiles in zip(
                    reactant_pair_smiles_tuples,
                    reaction_name_tuples,
                    reaction_product_smiles_tuples,
                ):
                    reaction_smarts = AllChem.ReactionFromSmarts(
                        "{}>>{}".format(
                            ".".join(reactant_pair_smiles), reaction_product_smiles
                        ),
                        useSmiles=True,
                    )
                    reaction_temperature = [
                        actionsession["actions"][0]["content"]["temperature"]["value"]
                        for actionsession in encoded_recipes[reaction_name]["recipes"][
                            "standard"
                        ]["actionsessions"]
                        if actionsession["type"] == "stir"
                    ][0]

                    reaction_id = createReactionModel(
                        method_id=method_id,
                        reaction_class=reaction_name,
                        reaction_number=1,
                        intramolecular=False,
                        recipe_type="standard",
                        reaction_temperature=reaction_temperature,
                        reaction_smarts=reaction_smarts,
                    )

                    createProductModel(
                        reaction_id=reaction_id,
                        product_smiles=reaction_product_smiles,
                    )

                    for reactant_smi in reactant_pair_smiles:
                        previousreactionqueryset = checkPreviousReactionProducts(
                            reaction_id=reaction_id,
                            smiles=reactant_smi,
                        )
                        if previousreactionqueryset:
                            reactant_id = createReactantModel(
                                reaction_id=reaction_id,
                                reactant_smiles=reactant_smi,
                                previous_reaction_product=True,
                            )
                            createCatalogEntryModel(
                                reactant_id=reactant_id,
                                previous_reaction_product=True,
                            )

                        if not previousreactionqueryset:
                            reactant_id = createReactantModel(
                                reaction_id=reaction_id,
                                reactant_smiles=reactant_smi,
                                previous_reaction_product=False,
                            )
                            catalog_entries = getExactSearch(smiles=reactant_smi)
                            if "results" in catalog_entries:
                                for catalog_entry in catalog_entries["results"]:
                                    createCatalogEntryModel(
                                        catalog_entry=catalog_entry,
                                        reactant_id=reactant_id,
                                        previous_reaction_product=False,
                                    )

    delete_tmp_file(csv_fp)

    return validate_dict, validated, project_info


@shared_task
def createOTScript(batchids: list, protocol_name: str):
    """ "
    Create otscripts and starting plates for a list of batch ids
    """
    task_summary = {}
    otprojectobj = OTProject()
    projectobj = Batch.objects.get(id=batchids[0]).project_id
    otprojectobj.project_id = projectobj
    otprojectobj.name = protocol_name
    otprojectobj.save()

    for batchid in batchids:
        reactionqueryset = getBatchReactions(batchid=batchid)
        # if reactionqueryset:
        #     allreactionobjs = [
        #         reactionobj
        #         for sublist in allreactionquerysets
        #         for reactionobj in sublist
        #     ]
        # reaction_ids = [reactionobj.id for reactionobj in reactionqueryset]
        actionsessionqueryset = getActionSessionQuerySet(reaction_ids=reactionqueryset)
        if not actionsessionqueryset:
            for reactionobj in reactionqueryset:
                reaction_id = reactionobj.id
                reactionclass = reactionobj.reactionclass
                target_id = reactionobj.method_id.target_id.id
                reactant_pair_smiles = reactionobj.reactants.all().values_list(
                    "smiles", flat=True
                )
                std_recipe = encoded_recipes[reactionclass]["recipes"][
                    "standard"
                ]  # NB need to include multiple recipes
                intramolecular_possible = encoded_recipes[reactionclass][
                    "intramolecular"
                ]

                if intramolecular_possible and len(reactant_pair_smiles) == 1:
                    intramolecular_possible = True
                else:
                    intramolecular_possible = False

                actionsessions = std_recipe["actionsessions"]

                CreateEncodedActionModels(
                    intramolecular=intramolecular_possible,
                    actionsessions=actionsessions,
                    target_id=target_id,
                    reaction_id=reaction_id,
                    reactant_pair_smiles=list(reactant_pair_smiles),
                    reaction_name=reactionclass,
                )

        batchtag = getBatchTag(batchid=batchid)
        otbatchprotocolqueryset = getOTBatchProtocolQuerySet(batch_id=batchid)

        if otbatchprotocolqueryset and otbatchprotocolqueryset[0].zipfile:
            otbatchprotocolobj = otbatchprotocolqueryset[0]
            otbatchprotocolobj.otproject_id = otprojectobj
            otbatchprotocolobj.celery_taskid = current_task.request.id
            otbatchprotocolobj.save()
            task_summary[batchid] = True
        else:
            otbatchprotocolobj = OTBatchProtocol()
            otbatchprotocolobj.batch_id = Batch.objects.get(id=batchid)
            otbatchprotocolobj.otproject_id = otprojectobj
            otbatchprotocolobj.celery_taskid = current_task.request.id
            otbatchprotocolobj.save()
            # maxsteps = findmaxlist(reactionqueryset=reactionqueryset)
            maxreactionnumber = getMaxReactionNumber(reactionqueryset=reactionqueryset)
            groupedreactionquerysets = groupReactions(
                reactionqueryset=reactionqueryset, maxreactionnumber=maxreactionnumber
            )
            for index, groupreactionqueryset in enumerate(groupedreactionquerysets):
                if index == 0:
                    reaction_ids = [reaction.id for reaction in groupreactionqueryset]
                    actionsessionqueryset = getActionSessionQuerySet(
                        reaction_ids=reaction_ids
                    )
                    sessionnumbers = getActionSessionSequenceNumbers(
                        actionsessionqueryset=actionsessionqueryset
                    )
                    groupedactionsessionsequences = getGroupedActionSessionSequences(
                        sessionnumbers=sessionnumbers,
                        actionsessionqueryset=actionsessionqueryset,
                    )
                    for groupactionsession in groupedactionsessionsequences:
                        actionsessiontypes = getActionSessionTypes(
                            actionsessionqueryset=groupactionsession
                        )
                        groupedactionsessiontypes = getGroupedActionSessionTypes(
                            actionsessiontypes=actionsessiontypes,
                            actionsessionqueryset=groupactionsession,
                        )
                        for groupactionsessiontype in groupedactionsessiontypes:
                            human_actionsessionqueryset = groupactionsessiontype.filter(
                                driver="human"
                            )
                            robot_actionsessionqueryset = groupactionsessiontype.filter(
                                driver="robot"
                            )
                            if human_actionsessionqueryset:
                                pass
                            if robot_actionsessionqueryset:

                                actionsession_ids = (
                                    robot_actionsessionqueryset.values_list(
                                        "id", flat=True
                                    )
                                )
                                session = CreateOTSession(
                                    reactionstep=index + 1,
                                    otbatchprotocolobj=otbatchprotocolobj,
                                    actionsessionqueryset=robot_actionsessionqueryset,
                                    groupreactionqueryset=groupreactionqueryset,
                                )

                                OTWrite(
                                    protocolname=batchtag,
                                    otsessionobj=session.otsessionobj,
                                    reaction_ids=reaction_ids,
                                    actionsession_ids=actionsession_ids,
                                )

                if index > 0:
                    groupreactiontodoqueryset = getReactionsToDo(
                        groupreactionqueryset=groupreactionqueryset
                    )
                    if len(groupreactiontodoqueryset) == 0:
                        break
                    else:
                        reaction_ids = [
                            reaction.id for reaction in groupreactiontodoqueryset
                        ]
                        actionsessionqueryset = getActionSessionQuerySet(
                            reaction_ids=reaction_ids
                        )
                        sessionnumbers = getActionSessionSequenceNumbers(
                            actionsessionqueryset=actionsessionqueryset
                        )
                        groupedactionsessionsequences = (
                            getGroupedActionSessionSequences(
                                sessionnumbers=sessionnumbers,
                                actionsessionqueryset=actionsessionqueryset,
                            )
                        )
                        for groupactionsession in groupedactionsessionsequences:
                            actionsessiontypes = getActionSessionTypes(
                                actionsessionqueryset=groupactionsession
                            )
                            groupedactionsessiontypes = getGroupedActionSessionTypes(
                                actionsessiontypes=actionsessiontypes,
                                actionsessionqueryset=groupactionsession,
                            )
                            for groupactionsessiontype in groupedactionsessiontypes:
                                human_actionsessionqueryset = (
                                    groupactionsessiontype.filter(driver="human")
                                )
                                robot_actionsessionqueryset = (
                                    groupactionsessiontype.filter(driver="robot")
                                )
                                if human_actionsessionqueryset:
                                    pass
                                if robot_actionsessionqueryset:
                                    actionsession_ids = (
                                        robot_actionsessionqueryset.values_list(
                                            "id", flat=True
                                        )
                                    )
                                    session = CreateOTSession(
                                        reactionstep=index + 1,
                                        otbatchprotocolobj=otbatchprotocolobj,
                                        actionsessionqueryset=robot_actionsessionqueryset,
                                        groupreactionqueryset=groupreactiontodoqueryset,
                                    )

                                    OTWrite(
                                        protocolname=batchtag,
                                        otsessionobj=session.otsessionobj,
                                        reaction_ids=reaction_ids,
                                        actionsession_ids=actionsession_ids,
                                    )

            createZipOTBatchProtocol = ZipOTBatchProtocol(
                otbatchprotocolobj=otbatchprotocolobj, batchtag=batchtag
            )

            if createZipOTBatchProtocol.errors:
                task_summary[batchid] = False
            else:
                task_summary[batchid] = True
    else:
        task_summary[batchid] = False

    return task_summary, otprojectobj.id


def getOTBatchProtocolQuerySet(batch_id: int) -> QuerySet[OTBatchProtocol]:
    """Gets the OT batch protocol queryset

    Parameters
    ----------
    batch_id: int
        The batch id to saerch for an associated OT protocol

    Returns
    -------
        The OT Batch protocol queryset
    """
    otbatchprotocolqueryset = OTBatchProtocol.objects.filter(batch_id=batch_id)
    return otbatchprotocolqueryset


def getActionSessionSequenceNumbers(
    actionsessionqueryset: QuerySet[ActionSession],
) -> list:
    """Set of action session sequence numbers

    Returns
    ------
    sessionnumbers: list
        The set of session numbers in an action session
        queryset eg. [1,2,3,4....n]
    """
    maxsessionnumber = actionsessionqueryset.aggregate(Max("sessionnumber"))[
        "sessionnumber__max"
    ]
    sessionnumbers = list(range(1, maxsessionnumber + 1))
    return sessionnumbers


def getActionSessionTypes(actionsessionqueryset: QuerySet[ActionSession]) -> QuerySet:
    """Set of action session types

    Returns
    ------
    actionsessiontypes: QuerySet
        The set of action session types in a queryset
        eg. ["reaction", "workup", "stir"]
    """
    actionsessiontypes = set(list(actionsessionqueryset.values_list("type", flat=True)))
    return actionsessiontypes


def getGroupedActionSessionSequences(
    sessionnumbers: list, actionsessionqueryset: QuerySet[ActionSession]
) -> list:
    """Group action sessions by sequence number

    Parameters
    ----------
    sessionnumbers: list
        The list of action session sequence numbers
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset to group by sequence number

    Returns
    -------
    groupedactionsessionsequences: list
        List of sublists of action sessions grouped by sequence number
    """
    groupedactionsessionsequences = []
    for sessionnumber in sessionnumbers:
        actionsessiongroup = actionsessionqueryset.filter(
            sessionnumber=sessionnumber
        ).order_by("-pk")
        groupedactionsessionsequences.append(actionsessiongroup)
    return groupedactionsessionsequences


def getGroupedActionSessionTypes(
    actionsessiontypes: QuerySet, actionsessionqueryset: QuerySet[ActionSession]
) -> list:
    """Group action sessions by type

    Parameters
    ----------
    actionsessiontypes: QuerySet
        The list of action session sequence numbers
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset to group by sequence number

    Returns
    -------
    groupedactionsessionquerysettypes: list
        List of sublists of action sessions grouped by types
    """
    groupedactionsessiontypes = []
    for actionsessiontype in actionsessiontypes:
        actionsessiongrouptype = actionsessionqueryset.filter(
            type=actionsessiontype
        ).order_by("-pk")
        if actionsessiongrouptype:
            groupedactionsessiontypes.append(actionsessiongrouptype)
    return groupedactionsessiontypes


def getActionSessionQuerySet(
    reaction_ids: QuerySet[Reaction],
    driver: str = None,
) -> QuerySet[ActionSession]:
    """Returns the action session wueryset for a type of driver
       (human or robot)

    Parameters
    ----------
    reactions_ids: QuerySet[Reaction]
        The reactions that the action session will excecute
    driver: str
        The optional main driver of the action session

    Returns
    -------
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset for a given driver
    """
    if driver:
        criterion1 = Q(reaction_id__in=reaction_ids)
        criterion2 = Q(driver=driver)
        actionsessionqueryset = ActionSession.objects.filter(criterion1 & criterion2)
        if actionsessionqueryset:
            return actionsessionqueryset
    if not driver:
        criterion1 = Q(reaction_id__in=reaction_ids)
        actionsessionqueryset = ActionSession.objects.filter(criterion1)
        if actionsessionqueryset:
            return actionsessionqueryset


def getPreviousObjEntries(queryset: list, obj: object):
    """Finds all previous objects relative to obj of queryset"""
    previousqueryset = queryset.filter(pk__lt=obj.pk).order_by("-pk")
    return previousqueryset


def checkPreviousReactionFailures(reactionobj: object):
    """Check if any previous reaction failures for a method"""
    reactionqueryset = getReactions(methodid=reactionobj.method_id.id)
    previousreactionqueryset = getPreviousObjEntries(
        queryset=reactionqueryset, obj=reactionobj
    )
    failedreactions = previousreactionqueryset.filter(success=False)
    if failedreactions.exists():
        return True
    else:
        return False


def checkNoMethodSteps(reactionobj):
    """Check no reaction steps in method is > 1"""
    methodobj = reactionobj.method_id
    noreactionsteps = methodobj.nosteps
    if noreactionsteps > 1:
        return True
    else:
        return False


def getReactionsToDo(groupreactionqueryset: QuerySet[Reaction]) -> QuerySet[Reaction]:
    """Get reactions that need to be done. Exclude those in methods that had
    failed previous reaction step

    Parameters
    ---------
    groupreactionqueryset: QuerySet[Reaction]
        The group of reactions to find to do based on if the previous reaction was successful

    Returns
    -------
    groupreactiontodoqueryset: QuerySet[Reaction]
        The reactions that need to be done
    """
    reactionstodo = []
    for reactionobj in groupreactionqueryset:
        if checkNoMethodSteps(reactionobj=reactionobj):
            if not checkPreviousReactionFailures(reactionobj=reactionobj):
                reactionstodo.append(reactionobj)
    groupreactiontodoqueryset = groupreactionqueryset.filter(
        reaction_id__in=reactionstodo
    )
    return groupreactiontodoqueryset


def getTargets(batch_ids: QuerySet[Batch]) -> QuerySet[Target]:
    targetqueryset = Target.objects.filter(batch_id__in=batch_ids).order_by("id")
    return targetqueryset


def getMethods(target_ids: QuerySet[Target]) -> QuerySet[Method]:
    methodqueryset = (
        Method.objects.filter(target_id__in=target_ids)
        .filter(otchem=True)
        .order_by("id")
    )
    return methodqueryset


def getReactions(method_ids: QuerySet[Method]) -> QuerySet[Reaction]:
    reactionqueryset = Reaction.objects.filter(method_id__in=method_ids).order_by("id")
    return reactionqueryset


def getBatchTag(batchid):
    batch_obj = Batch.objects.get(id=batchid)
    batchtag = batch_obj.batchtag
    return batchtag


def getBatchReactions(batchid: int) -> QuerySet[Reaction]:
    targetqueryset = getTargets(batch_ids=[batchid])
    if targetqueryset:
        methodqueryset = getMethods(target_ids=targetqueryset)
        if methodqueryset:
            reactionqueryset = getReactions(method_ids=methodqueryset)
            if reactionqueryset:
                return reactionqueryset


def getMaxReactionNumber(reactionqueryset: QuerySet[Reaction]) -> int:
    """Get the maximum number of reaction steps in a reaction queryset

    Parameters
    ----------
    reactionqueryset: QuerySet[Reaction]
        The reaction queryset to get the max number of reaction steps for

    Returns
    -------
    maxreactionnumber: int
        The maximum reaction number in a set of reactions

    """

    maxreactionnumber = reactionqueryset.aggregate(Max("number"))["number__max"]
    return maxreactionnumber


def groupReactions(reactionqueryset: QuerySet[Reaction], maxreactionnumber: int):
    """
    Groups reactionqueries into first reactions, second reactions and so on
    """

    groupedreactionquerysets = []
    for i in range(1, maxreactionnumber + 1):
        reactionnumberqueryset = (
            reactionqueryset.filter(number=i).distinct().order_by("id")
        )
        if reactionnumberqueryset:
            groupedreactionquerysets.append(reactionnumberqueryset)
    return groupedreactionquerysets


class ZipOTBatchProtocol(object):
    """
    Creates a ZipBatchProtocol object for writing compound order csvs and otscripts
    for a batch
    """

    def __init__(self, otbatchprotocolobj: object, batchtag: str):
        """
        zipOTBatchProtocol constructor
        Args:
            otbatchprotocolobj (Django object): OT batch protocol object created for a batch
            batchtag (str): Batch tag for naming zip file
        """
        self.errors = {"function": [], "errorwarning": []}
        self.otbatchprotocolobj = otbatchprotocolobj
        self.mediaroot = settings.MEDIA_ROOT
        self.otsessionqueryset = self.getOTSessionQuerySet()
        self.zipfn = "batch-{}-protocol.zip".format(batchtag)
        self.ziptmpfp = os.path.join(settings.MEDIA_ROOT, "tmp", "batchprotocoltmp.zip")
        self.ziparchive = ZipFile(self.ziptmpfp, "w")

        for otsession_obj in self.otsessionqueryset:
            solventprepqueryset = self.getSolventPrepQuerySet(
                otsessionobj=otsession_obj
            )
            compoundorderqueryset = self.getCompoundOrderQuerySet(
                otsessionobj=otsession_obj
            )
            otscriptqueryset = self.getOTScriptQuerySet(otsessionobj=otsession_obj)

            if solventprepqueryset:
                for solventprepobj in solventprepqueryset:
                    filepath = self.getSolventPrepFilePath(
                        solventprepobj=solventprepobj
                    )
                    destdir = "solventprep"
                    self.writeZip(destdir=destdir, filepath=filepath)

            if compoundorderqueryset:
                for compoundorderobj in compoundorderqueryset:
                    filepath = self.getCompoundOrderFilePath(
                        compoundorderobj=compoundorderobj
                    )
                    destdir = "compoundorders"
                    self.writeZip(destdir=destdir, filepath=filepath)

            if otscriptqueryset:
                for otscriptobj in otscriptqueryset:
                    filepath = self.getOTScriptFilePath(otscriptobj=otscriptobj)
                    destdir = "otscripts"
                    self.writeZip(destdir=destdir, filepath=filepath)

        self.ziparchive.close()
        self.writeZipToMedia()
        self.deleteTmpZip()

    def addWarning(self, function: str, errorwarning: str):
        self.errors["function"].append(function)
        self.errors["errorwarning"].append(errorwarning)

    def getOTSessionQuerySet(self):
        """Retrieve OTSession model queryset"""
        otsessionqueryset = OTSession.objects.filter(
            otbatchprotocol_id=self.otbatchprotocolobj
        )

        if not otsessionqueryset:
            self.addWarning(
                function=self.getOTSessionQuerySet.__name__,
                errorwarning="No queryset found",
            )
        else:
            return otsessionqueryset

    def getSolventPrepQuerySet(self, otsessionobj: object):
        """Retrieve SolventPrep model queryset
        Args:
            otsessionobj (Django obj): OTSession Django object
        """
        solventprepqueryset = SolventPrep.objects.filter(otsession_id=otsessionobj)

        if not solventprepqueryset:
            self.addWarning(
                function=self.getSolventPrepQuerySet.__name__,
                errorwarning="No queryset found",
            )
            return None
        else:
            return solventprepqueryset

    def getCompoundOrderQuerySet(self, otsessionobj: object):
        """Retrieve CompoundOrder model queryset
        Args:
            otsessionobj (Django obj): OTSession Django object
        """
        compoundorderqueryset = CompoundOrder.objects.filter(otsession_id=otsessionobj)

        if not compoundorderqueryset:
            self.addWarning(
                function=self.getCompoundOrderQuerySet.__name__,
                errorwarning="No queryset found",
            )
        else:
            return compoundorderqueryset

    def getOTScriptQuerySet(self, otsessionobj: object):
        """Retrieve OTScript model queryset
        Args:
            otsessionobj (Django obj): OTSession Django object
        """
        otscriptqueryset = OTScript.objects.filter(otsession_id=otsessionobj)

        if not otscriptqueryset:
            self.addWarning(
                function=self.getOTScriptQuerySet.__name__,
                errorwarning="No queryset found",
            )
        else:
            return otscriptqueryset

    def getSolventPrepFilePath(self, solventprepobj: object):
        """Retrieve SolventPrep csv file path
        Args:
            solventprepobj (Django obj): SolventPrep Django object
        """
        filepath = os.path.join(self.mediaroot, solventprepobj.solventprepcsv.name)
        return filepath

    def getCompoundOrderFilePath(self, compoundorderobj: object):
        """Retrieve CompoundOrder csv file path
        Args:
            compoundorderobj (Django obj): OTSession Django object
        """
        filepath = os.path.join(self.mediaroot, compoundorderobj.ordercsv.name)
        return filepath

    def getOTScriptFilePath(self, otscriptobj: object):
        """Retrieve OTScript Python file path
        Args:
            otscriptobj (Django obj): OTScript Django object
        """
        filepath = os.path.join(self.mediaroot, otscriptobj.otscript.name)
        return filepath

    def writeZip(self, destdir: str, filepath: str):
        """Add the requested file to the zip archive.

        Args:
            destdir (str): directory to write file to in ziparchive
            filepath (str): filepath from record
        """
        arcname = os.path.join(destdir, filepath.split("/")[-1])
        self.ziparchive.write(filename=filepath, arcname=arcname)

    def writeZipToMedia(self):
        """Write the ziparchive to medida for the
        OTBatchProtocol Django object
        """
        zf = open(self.ziptmpfp, "rb")
        self.otbatchprotocolobj.zipfile.save(self.zipfn, zf)
        self.otbatchprotocolobj.save()

    def deleteTmpZip(self):
        """ "Delete the temporary zip archive created by ZipFile"""
        os.remove(self.ziptmpfp)


@shared_task
def canonicalizeSmiles(csvfile: str = None, smiles: list = None):
    """ "
    Canonicalizes smiles from csv file uploaded from frontend
    """
    validated = True

    if csvfile:
        csvdf = pd.read_csv(csvfile, encoding="utf8")
        smiles = [smi for smi in csvdf["SMILES"]]
        delete_tmp_file(csvfile)
    if smiles:
        smiles = smiles

    molcheck = [Chem.MolFromSmiles(smi) for smi in smiles]

    if None not in molcheck:
        canonicalizedsmiles = [canonSmiles(smi) for smi in smiles]
        return validated, canonicalizedsmiles
    else:
        validated = False
        indexerrors = [i for i, v in enumerate(molcheck) if v is None]
        errorsummary = "There was an error with the smiles csv at index: {}".format(
            indexerrors
        )
        return validated, errorsummary
