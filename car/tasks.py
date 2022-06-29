"""CAR's celery tasks"""
from __future__ import annotations
from celery import shared_task, current_task
from django.conf import settings
from django.db.models import QuerySet

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

from .manifold.apicalls import getManifoldRetrosynthesis
from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import getAddtionOrder, canonSmiles

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

    return (validate_dict, validated, project_info, csv_fp, uploaded_dict)


@shared_task
def uploadManifoldReaction(validate_output):

    validate_dict, validated, project_info, csv_fp, uploaded_dict = validate_output
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

            for smiles, mass in zip(group["targets"], group["amount-required-mg"]):

                target_id = createTargetModel(
                    batch_id=batch_id,
                    smiles=smiles,
                    mass=mass,
                )

                retrosynthesis_result = getManifoldRetrosynthesis(smiles=smiles)
                routes = retrosynthesis_result["routes"]

                if not routes:
                    continue
                else:
                    first_route = routes[0]

                    # Check if target is in a catalogue and create catalog entries if it is
                    if first_route["molecules"][0]["isBuildingBlock"]:
                        catalog_entries = first_route["molecules"][0]["catalogEntries"]
                        for catalog_entry in catalog_entries:
                            createCatalogEntryModel(
                                catalog_entry=catalog_entry, target_id=target_id
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
                                target_id=target_id, nosteps=no_steps, otchem=False
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
                                    product_smiles=product_smiles,
                                )

                                for reactant_smi in reactant_smiles:
                                    reactant_id = createReactantModel(
                                        reaction_id=reaction_id,
                                        reactant_smiles=reactant_smi,
                                    )
                                    catalog_entries = [
                                        molecule["catalogEntries"]
                                        for molecule in route["molecules"]
                                        if molecule["smiles"] == reactant_smi
                                    ][0]
                                    for catalog_entry in catalog_entries:
                                        createCatalogEntryModel(
                                            catalog_entry=catalog_entry,
                                            reactant_id=reactant_id,
                                        )

                        if len(encoded_reactions_found) == no_steps:
                            method_id = createMethodModel(
                                target_id=target_id, nosteps=no_steps, otchem=True
                            )

                            for reaction in reversed(reactions):
                                reaction_name = reaction["name"]
                                reactant_smiles = reaction["reactantSmiles"]
                                product_smiles = reaction["productSmiles"]

                                std_recipe = encoded_recipes[reaction_name]["recipes"][
                                    "standard"
                                ]  # NB need to include multiple recipes
                                recipe_rxn_smarts = encoded_recipes[reaction_name][
                                    "reactionSMARTS"
                                ]
                                intramolecular_possible = encoded_recipes[
                                    reaction_name
                                ]["intramolecular"]

                                if (
                                    len(reactant_smiles) == 1
                                    and intramolecular_possible
                                ):
                                    intramolecular = True
                                    actionsessions = std_recipe["actionsessions"]
                                    stir_action = std_recipe["actionsessions"]["stir"][
                                        "actions"
                                    ][0]
                                    reaction_temperature = stir_action["content"][
                                        "temperature"
                                    ]["value"]
                                    reactant_smiles_ordered = reactant_smiles
                                else:
                                    intramolecular = False
                                    actionsessions = std_recipe["actionsessions"]
                                    stir_action = std_recipe["actionsessions"]["stir"][
                                        "actions"
                                    ][0]
                                    reaction_temperature = stir_action["content"][
                                        "temperature"
                                    ]["value"]
                                    reactant_smiles_ordered = getAddtionOrder(
                                        product_smi=product_smiles,
                                        reactant_SMILES=reactant_smiles,
                                        reaction_SMARTS=recipe_rxn_smarts,
                                    )
                                    if not reactant_smiles_ordered:
                                        continue

                                reaction_smarts = AllChem.ReactionFromSmarts(
                                    "{}>>{}".format(
                                        ".".join(reactant_smiles_ordered),
                                        product_smiles,
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
                                    product_smiles=product_smiles,
                                )

                                for reactant_smi in reactant_smiles_ordered:
                                    reactant_id = createReactantModel(
                                        reaction_id=reaction_id,
                                        reactant_smiles=reactant_smi,
                                    )

                                    catalog_entries = [
                                        molecule["catalogEntries"]
                                        for molecule in route["molecules"]
                                        if molecule["smiles"] == reactant_smi
                                    ][0]
                                    # Get lowest cost/leadtime/preferred vendor - TO DO!!!!
                                    for catalog_entry in catalog_entries:
                                        createCatalogEntryModel(
                                            catalog_entry=catalog_entry,
                                            reactant_id=reactant_id,
                                        )

                                CreateEncodedActionModels(
                                    intramolecular=intramolecular,
                                    actionsessions=actionsessions,
                                    target_id=target_id,
                                    reaction_id=reaction_id,
                                    reactant_pair_smiles=reactant_smiles_ordered,
                                    reaction_name=reaction_name,
                                )

        delete_tmp_file(csv_fp)

        return validate_dict, validated, project_info


@shared_task
def uploadCustomReaction(validate_output):

    validate_dict, validated, project_info, csv_fp, uploaded_dict = validate_output
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

            for reactant_pair_smiles, reaction_name, smiles, mass in zip(
                group["reactant-pair-smiles"],
                group["reaction-name"],
                group["target-smiles"],
                group["amount-required-mg"],
            ):
                reaction_smarts = AllChem.ReactionFromSmarts(
                    "{}>>{}".format(".".join(reactant_pair_smiles), smiles),
                    useSmiles=True,
                )

                target_id = createTargetModel(
                    batch_id=batch_id,
                    smiles=smiles,
                    mass=mass,
                )

                recipes = encoded_recipes[reaction_name]["recipes"]

                actions = recipes["Standard"]["actions"]
                stir_action = [
                    action for action in actions if action["type"] == "stir"
                ][0]
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
                    product_smiles=smiles,
                )

                CreateEncodedActionModels(
                    actions=actions,
                    target_id=target_id,
                    reaction_id=reaction_id,
                    reactant_pair_smiles=reactant_pair_smiles,
                    reaction_name=reaction_name,
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
        allreactionquerysets = getBatchReactions(batchid=batchid)
        if allreactionquerysets:
            otbatchprotocolobj = OTBatchProtocol()
            otbatchprotocolobj.batch_id = Batch.objects.get(id=batchid)
            otbatchprotocolobj.otproject_id = otprojectobj
            otbatchprotocolobj.celery_taskid = current_task.request.id
            otbatchprotocolobj.save()
            batchtag = getBatchTag(batchid=batchid)
            maxsteps = findmaxlist(allreactionquerysets=allreactionquerysets)
            groupedreactionquerysets = groupReactions(
                allreactionquerysets=allreactionquerysets, maxsteps=maxsteps
            )
            for index, reactiongroup in enumerate(groupedreactionquerysets):
                if index == 0:
                    reaction_ids = [reaction.id for reaction in reactiongroup]
                    human_actionsessionqueryset = getActionSessions(
                        driver="human", reaction_ids=reaction_ids
                    )
                    robot_actionsessionqueryset = getActionSessions(
                        driver="robot", reaction_ids=reaction_ids
                    )
                    if human_actionsessionqueryset:
                        pass
                    if robot_actionsessionqueryset:
                        actionsessiontypes = getActionSessionTypes()
                        for actionsessiontype in actionsessiontypes:
                            type_actionsessionqueryset = getActionSessionByType(
                                actionsessiontype=actionsessiontype,
                                actionsessionqueryset=robot_actionsessionqueryset,
                            )
                            actionsession_ids = type_actionsessionqueryset.values_list(
                                "id", flat=True
                            )

                            session = CreateOTSession(
                                reactionstep=index + 1,
                                otbatchprotocolobj=otbatchprotocolobj,
                                actionsessionqueryset=type_actionsessionqueryset,
                                reactiongrouplist=reactiongroup,
                            )

                            OTWrite(
                                protocolname=batchtag,
                                otsessionobj=session.otsessionobj,
                                reaction_ids=reaction_ids,
                                actionsession_ids=actionsession_ids,
                            )

                if index > 0:
                    reactiongrouptodo = getReactionsToDo(reactiongroup=reactiongroup)
                    if len(reactiongrouptodo) == 0:
                        break
                    else:
                        reaction_ids = [reaction.id for reaction in reactiongrouptodo]
                        human_actionsessionqueryset = getActionSessions(
                            driver="human", reaction_ids=reaction_ids
                        )
                        robot_actionsessionqueryset = getActionSessions(
                            driver="robot", reaction_ids=reaction_ids
                        )
                        if human_actionsessionqueryset:
                            pass
                        if robot_actionsessionqueryset:
                            actionsessiontypes = getActionSessionTypes()
                            for actionsessiontype in actionsessiontypes:
                                type_actionsessionqueryset = getActionSessionByType(
                                    actionsessiontype=actionsessiontype,
                                    actionsessionqueryset=robot_actionsessionqueryset,
                                )
                                actionsession_ids = (
                                    type_actionsessionqueryset.values_list(
                                        "id", flat=True
                                    )
                                )

                                session = CreateOTSession(
                                    reactionstep=index + 1,
                                    otbatchprotocolobj=otbatchprotocolobj,
                                    actionsessionqueryset=type_actionsessionqueryset,
                                    reactiongrouplist=reactiongrouptodo,
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


def getActionSessions(driver: str, reaction_ids: list[int]) -> QuerySet[ActionSession]:
    """Returns the action session wueryset for a type of driver
       (human or robot)

    Parameters
    ----------
    driver: str
        The main driver of the action session
    reactions_ids: list[int]
        The reactions that the action session will excecute

    Returns
    -------
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset for a given driver
    """
    actionsessionqueryset = ActionSession.objects.filter(
        reaction_id__in=reaction_ids, driver=driver
    )
    return actionsessionqueryset


def getActionSessionByType(
    actionsessiontype: str, actionsessionqueryset: QuerySet[ActionSession]
) -> QuerySet[ActionSession]:
    """Returns the action session queryset for a type of
       action session eg. reaction, stir, analyse etc

    Parameters
    ----------
    actionsessiontype: str
        The actionsession type to filter the queryset
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset to filter for a given action session type

    Returns
    -------
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset for a given action session type
    """
    actionsessionqueryset = actionsessionqueryset.filter(type=actionsessiontype)
    return actionsessionqueryset


def getActionSessionTypes() -> list:
    """Returns the action session types eg. reaction, stir, analyse

    Returns
    -------
    actionsessiontypes: list
        The action session types
    """
    actionsessiontypes = [choice[0] for choice in ActionSession.type.field.choices]
    return actionsessiontypes


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
    if failedreactions:
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


def getReactionsToDo(reactiongroup):
    """get reactions that need to be done. Exclude those in methods that had
    failed previous reaction step
    """
    reactiongrouptodo = []
    for reactionobj in reactiongroup:
        if checkNoMethodSteps(reactionobj=reactionobj):
            if not checkPreviousReactionFailures(reactionobj=reactionobj):
                reactiongrouptodo.append(reactionobj)
    return reactiongrouptodo


def getTargets(batchid):
    targetqueryset = Target.objects.filter(batch_id=batchid).order_by("id")
    return targetqueryset


def getMethods(targetid):
    methodqueryset = (
        Method.objects.filter(target_id=targetid).filter(otchem=True).order_by("id")
    )
    return methodqueryset


def getReactions(methodid):
    reactionqueryset = Reaction.objects.filter(method_id=methodid).order_by("id")
    return reactionqueryset


def getBatchTag(batchid):
    batch_obj = Batch.objects.get(id=batchid)
    batchtag = batch_obj.batchtag
    return batchtag


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
        self.otbatchprotocolobj = otbatchprotocolobj
        self.mediaroot = settings.MEDIA_ROOT
        self.otsessionqueryset = self.getOTSessionQuerySet()
        self.zipfn = "batch-{}-protocol.zip".format(batchtag)
        self.ziptmpfp = os.path.join(settings.MEDIA_ROOT, "tmp", "batchprotocoltmp.zip")
        self.ziparchive = ZipFile(self.ziptmpfp, "w")
        self.errors = {"function": [], "errorwarning": []}

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
