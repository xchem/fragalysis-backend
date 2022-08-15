# this is vunrable to python injection by the lack of checking of metadata inputs
# this opens and closes files frequently, could be improved by creating string to hold the file data before writing to file once

#   NOTE: in this file "humanread" referes to the comments above each line/set of lines of ot code, human readable is a list of all the comments in format [oporator (human/ot), comment]
"""Create otWrite object"""
from __future__ import annotations
from django.core.files.storage import default_storage
from django.conf import settings
from django.db.models import QuerySet, Q
import os
from graphene_django import DjangoObjectType

from car.recipebuilder.encodedrecipes import encoded_recipes
from car.models import (
    ActionSession,
    AddAction,
    Column,
    ExtractAction,
    MixAction,
    OTSession,
    Reaction,
)

from car.models import (
    Product,
    Pipette,
    TipRack,
    Plate,
    Well,
    OTScript,
)
import math
import inspect
import logging

from .labwareavailable import labware_plates

logger = logging.getLogger(__name__)


class OTWrite(object):
    """ "
    Creates a otWrite object for generating an OT protocol
    script
    """

    def __init__(
        self,
        protocolname: str,
        otsessionobj: OTSession,
        reaction_ids: list,
        actionsession_ids: list,
        apiLevel: str = "2.9",
    ):
        """Initiates an otWrite object

        Parameters
        ----------
        protocolname: str
            The name of the OT protcol
        otsessionobj: OTSession
            The OT session the protocol is related to
        apiLevel: str = "2.9"
            The OpenTrons API version used
        reaction_ids: list
            The reactions that the an OT script will be generated for
        actionsession_ids: list
            The action session ids that the protocol is being written for
        """
        self.reactionstep = otsessionobj.reactionstep
        self.otsessionobj = otsessionobj
        self.otsession_id = otsessionobj.id
        self.otsessiontype = otsessionobj.sessiontype
        self.protocolname = protocolname
        self.apiLevel = apiLevel
        self.reaction_ids = reaction_ids
        self.actionsession_ids = actionsession_ids

        self.tiprackqueryset = self.getTipRacks()
        self.platequeryset = self.getPlates()
        self.pipetteobj = self.getPipette()
        self.pipettename = self.pipetteobj.name
        self.filepath, self.filename = self.createFilePath()
        self.setupScript()
        self.setupTipRacks()
        self.setupPlates()
        self.setupPipettes()

        if self.otsessiontype == "reaction":
            self.writeReactionSession()
        if self.otsessiontype == "workup":
            self.writeWorkUpSession()
        if self.otsessiontype == "analyse":
            self.writeAnalyseSession()

    def writeReactionSession(self):
        actionsessionqueryset = self.getActionSessionQuerySet()
        self.writeReactionActions(actionsessionqueryset=actionsessionqueryset)
        self.createOTScriptModel()

    def writeWorkUpSession(self):
        actionsessionqueryset = self.getActionSessionQuerySet()
        self.writeWorkUpActions(actionsessionqueryset=actionsessionqueryset)
        self.createOTScriptModel()

    def writeAnalyseSession(self):
        actionsessionqueryset = self.getActionSessionQuerySet()
        self.writeAnalyseActions(actionsessionqueryset=actionsessionqueryset)
        self.createOTScriptModel()

    def getActionSessionQuerySet(self) -> QuerySet[ActionSession]:
        """Get action session queryset for reaction_ids and actionsession_ids

        Returns
        -------
        addactionqueryset: QuerySet[AddAction]
            The add actions related to the reaction and action sessions
        """
        criterion1 = Q(reaction_id__in=self.reaction_ids)
        criterion2 = Q(id__in=self.actionsession_ids)

        actionsessionqueryset = ActionSession.objects.filter(
            criterion1 & criterion2
        ).order_by("id")
        return actionsessionqueryset

    def getAddActionQuerySet(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
        actionnumber: int = None,
    ) -> QuerySet[AddAction]:
        """Get add actions queryset for reaction_id

        Parameters
        ----------
        reaction_ids: list
            The reactions to search for related add actions
        actionsession_ids: list
            Optional action session ids to match add actions with
        actionsessiontype: str
            The optional action session type to look for add actions for
        actionnumber: int
            The optional number the action is executed in a process

        Returns
        -------
        addactionqueryset: QuerySet[AddAction]
            The add actions related to the reaction
        """

        if actionsession_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__in=actionsession_ids)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return addactionqueryset
        if actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return addactionqueryset
        if actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            criterion3 = Q(number=actionnumber)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2 & criterion3
            ).order_by("id")
            return addactionqueryset

    def getExtractActionQuerySet(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
        actionnumber: int = None,
    ) -> QuerySet[ExtractAction]:
        """Get extract actions queryset for reaction_ids

        Parameters
        ----------
        reaction_ids: list
            The reactions to search for related extract actions
        actionsession_ids: list
            Optional action session ids to match extract actions with
        actionsessiontype: str
            The optional action session type to look for extract actions for
        actionnumber: int
            The optional number the action is executed in a process

        Returns
        -------
        extractactionqueryset: QuerySet[ExtractAction]
            The extract actions related to the reaction
        """

        if actionsession_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__in=actionsession_ids)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return extractactionqueryset
        if actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return extractactionqueryset
        if actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            criterion3 = Q(number=actionnumber)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1 & criterion2 & criterion3
            ).order_by("id")
            return extractactionqueryset

    def getMixActionQuerySet(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
        actionnumber: int = None,
    ) -> QuerySet[MixAction]:
        """Get mix actions queryset for reaction_ids

        Parameters
        ----------
        reaction_ids: list
            The reactions to search for related mix actions
        actionsession_ids: list
            Optional action session ids to match mix actions with
        actionsessiontype: str
            The optional action session type to look for extract actions for
        actionnumber: int
            The optional number the action is executed in a process

        Returns
        -------
        mixactionqueryset: QuerySet[MixAction]
            The mix actions related to the reaction
        """

        if actionsession_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__in=actionsession_ids)
            mixactionqueryset = MixAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return mixactionqueryset
        if actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            mixactionqueryset = MixAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return mixactionqueryset
        if actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            criterion3 = Q(number=actionnumber)
            mixactionqueryset = MixAction.objects.filter(
                criterion1 & criterion2 & criterion3
            ).order_by("id")
            return mixactionqueryset

    def getAddAction(self, reaction_id: int) -> Product:
        """Gets the product for a reaction

        Parameters
        ----------
        reaction_id: int
            The reaction objects id to search for a product

        Returns
        -------
        productobj: Product
            The reaction's product
        """
        productobj = Product.objects.get(reaction_id=reaction_id)
        return productobj

    def getProduct(self, reaction_id: int) -> Product:
        """Gets the product for a reaction

        Parameters
        ----------
        reaction_id: int
            The reaction objects id to search for a product

        Returns
        -------
        productobj: Product
            The reaction's product
        """
        productobj = Product.objects.get(reaction_id=reaction_id)
        return productobj

    def getReaction(self, reaction_id: int) -> Reaction:
        """Gets the reaction object

        Parameters
        ----------
        reaction_id: int
            The reaction objects id

        Returns
        -------
        reactionobj: Reaction
            The reaction object
        """
        reactionobj = Reaction.objects.get(id=reaction_id)
        return reactionobj

    def getPlates(self) -> QuerySet[Plate]:
        """Gets plates for an OT session

        Returns
        -------
        platequeryset: QuerySet[Plate]
            The plates linked to the OT session
        """
        platequeryset = Plate.objects.filter(otsession_id=self.otsession_id).order_by(
            "id"
        )
        return platequeryset

    def getPlateObj(self, plateid: int) -> Plate:
        """Gets the plate object

        Parameters
        ----------
        plateid: int
            The id of the plate to search for

        Returns
        -------
        plateobj: Plate
            The plate object
        """
        plateobj = Plate.objects.filter(id=plateid)[0]
        return plateobj

    def getTipRacks(self) -> QuerySet[TipRack]:
        """Get the tip racks for an OT Session

        Returns
        -------
        tiprackqueryset: QuerySet[Tiprack]
            The tipracks linked to an OT Session
        """
        tipracksqueryset = TipRack.objects.filter(
            otsession_id=self.otsession_id
        ).order_by("id")
        return tipracksqueryset

    def getPipette(self) -> Pipette:
        """Get the pipette for an OT session

        Returns
        -------
        pipetteobj: Pipette
            The pipette used for an OT session
        """
        pipetteobj = Pipette.objects.get(otsession_id=self.otsession_id)
        return pipetteobj

    def getProductSmiles(self, reaction_id: int) -> str:
        """Get the product smiles

        Parameters
        ----------
        reaction_id: int
            The reaction id linked to a product

        Returns
        -------
        productobj.smiles: str
            The SMILES of the product of  reaction
        """
        productobj = Product.objects.filter(reaction_id=reaction_id)[0]
        return productobj.smiles

    def getPreviousObjEntries(
        self, queryset: QuerySet, obj: DjangoObjectType
    ) -> QuerySet:
        """Finds all previous Django model object relative to the Django model
           object in a queryset

        Parameters
        ----------
        queryset: QuerySet
            The queryset to search for previous entries
        obj: DjangoObjectType
            The object that you want to find all previous object entries relative to

        Returns
        -------
        previousqueryset: QuerySet
            The previous Django model objects as a queryset
        """
        previousqueryset = queryset.filter(pk__lt=obj.pk).order_by("-pk")
        return previousqueryset

    def getReactionQuerySet(
        self, reaction_ids: list = None, method_id: int = None
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
            reactionqueryset = Reaction.objects.filter(id__in=reaction_ids)
        if method_id:
            reactionqueryset = Reaction.objects.filter(method_id=method_id)
        return reactionqueryset

    def getColumnQuerySet(
        self,
        columntype: str,
        reactionclass: str,
    ) -> QuerySet[Column]:
        """Get column queryset for column type and reactionclass

        Parameters
        ----------
        columntype: str
            The type of plate the column is on eg. reaction, workup1, workup2, analyse
        reactionclass: str
            The reaction class that occupy the columns

        Returns
        -------
        columnqueryset: QuerySet[Column]
            The columns related to the column type and reaction class
        """
        print(self.otsession_id, columntype, reactionclass)
        criterion1 = Q(otsession_id=self.otsession_id)
        criterion2 = Q(type=columntype)
        criterion3 = Q(reactionclass=reactionclass)
        columnqueryset = Column.objects.filter(
            criterion1 & criterion2 & criterion3
        ).order_by("id")
        return columnqueryset

    def checkPreviousReactionProduct(
        self, reaction_id: int, smiles: str
    ) -> list[Reaction]:
        """Checks if any previous reactions had a product matching the smiles

        Parameters
        ----------
        reaction_id: int
            The reaction id of the Django model object to search for
            all relative previous reactions objects. The previous reactions may
            have products that are this reaction's reactant input
        smiles: str
            The SMILES of the reaction's reactant and previous reaction products

        Returns
        -------
        previousproductmatches: list[Reactant]
            The list of reactions whose product matches the query SMILES
        """

        reactionobj = self.getReaction(reaction_id=reaction_id)
        reactionqueryset = self.getReactionQuerySet(method_id=reactionobj.method_id.id)
        prevreactionqueryset = self.getPreviousObjEntries(
            queryset=reactionqueryset, obj=reactionobj
        )
        previousproductmatches = []
        if prevreactionqueryset:
            for reactionobj in prevreactionqueryset:
                productobj = self.getProduct(reaction_id=reactionobj)
                if productobj.smiles == smiles:
                    previousproductmatches.append(reactionobj)
        return previousproductmatches

    def createFilePath(self):
        """Creates the OT protcol script filepath"""
        filename = (
            "{}-session-ot-script-batch-{}-reactionstep{}-sessionid-{}.txt".format(
                self.otsessiontype,
                self.protocolname,
                self.reactionstep,
                self.otsession_id,
            )
        )
        path = "tmp/" + filename
        filepath = str(os.path.join(settings.MEDIA_ROOT, path))
        return filepath, filename

    def createOTScriptModel(self):
        """Creates an OT Script model object"""
        script = open(self.filepath, "a")
        script.close()

        otscriptobj = OTScript()
        otscriptobj.otsession_id = self.otsessionobj
        otscriptfile = open(self.filepath, "rb")
        otscriptfn = default_storage.save(
            "otscripts/{}.py".format(self.filename.strip(".txt")), otscriptfile
        )
        otscriptobj.otscript = otscriptfn
        otscriptfile.close()
        otscriptobj.save()

    def findSolventPlateWellObj(self, solvent: str, transfervolume: float) -> list:
        """Finds solvent well for diluting a previous reaction steps product

        Parameters
        ----------
        solvent: str
            The solvent needed for diluting the contents of a well
        transfervolume: float
            The volume of solvent required for dilution

        Returns
        -------
        wellinfo: list
            The list of wells found along with volume to transfer from the
            well
        """
        wellinfo = []
        try:
            wellobjs = Well.objects.filter(
                otsession_id=self.otsession_id,
                solvent=solvent,
                available=True,
                type="solvent",
            ).order_by("id")
            for wellobj in wellobjs:
                areclose = self.checkVolumeClose(volume1=transfervolume, volume2=0.00)
                if areclose:
                    break
                wellvolumeavailable = self.getWellVolumeAvailable(wellobj=wellobj)
                if wellvolumeavailable > 0:
                    if wellvolumeavailable >= transfervolume:
                        self.updateWellVolume(
                            wellobj=wellobj, transfervolume=transfervolume
                        )
                        wellinfo.append([wellobj, transfervolume])
                        transfervolume = 0.00
                    if wellvolumeavailable < transfervolume:
                        self.updateWellVolume(
                            wellobj=wellobj, transfervolume=wellvolumeavailable
                        )
                        wellinfo.append([wellobj, wellvolumeavailable])
                        transfervolume = transfervolume - wellvolumeavailable
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            print(e)
            print(solvent)
        return wellinfo

    def findStartingPlateWellObj(
        self,
        reaction_id: int,
        smiles: str,
        solvent: str,
        concentration: float,
        transfervolume: float,
    ) -> list:
        """Finds starting plate well for executing an add action

        Parameters
        ----------
        reaction_id: int
            The reaction's id used in the transfer
        smiles: str
            The SMILES of the starting material needed to be tranferred
        solvent: str
            The solvent used to prepare the starting material
        concentration: float
            The concentration of the starting material
        transfervolume: float
            The volume of starting material needed for the transfer

        Returns
        -------
        wellinfo: list
            The list of wells found along with volume to transfer from the
            well
        """
        previousreactionobjs = self.checkPreviousReactionProduct(
            reaction_id=reaction_id, smiles=smiles
        )
        wellinfo = []
        if previousreactionobjs:
            criterion1 = Q(otsession_id=self.otsession_id)
            criterion2 = Q(reaction_id=previousreactionobjs[0])
            criterion3 = Q(smiles=smiles)
            criterion4 = Q(type="reaction")
            criterion5 = Q(type="workup1")
            criterion6 = Q(type="workup2")
            criterion7 = Q(type="workup3")
            try:
                wellobj = Well.objects.get(
                    criterion1
                    & criterion2
                    & criterion3
                    & (criterion5 | criterion6 | criterion7)
                )
                wellinfo.append([previousreactionobjs, wellobj, transfervolume])
            except:
                wellobj = Well.objects.get(
                    criterion1 & criterion2 & criterion3 & criterion4
                )
                wellinfo.append([previousreactionobjs, wellobj, transfervolume])
        else:
            try:
                wellobjects = Well.objects.filter(
                    otsession_id=self.otsession_id,
                    smiles=smiles,
                    solvent=solvent,
                    concentration=concentration,
                    available=True,
                    type="startingmaterial",
                ).order_by("id")
                for wellobj in wellobjects:
                    areclose = self.checkVolumeClose(
                        volume1=transfervolume, volume2=0.00
                    )
                    if areclose:
                        break
                    wellvolumeavailable = self.getWellVolumeAvailable(wellobj=wellobj)
                    if wellvolumeavailable > 0:
                        if wellvolumeavailable >= transfervolume:
                            self.updateWellVolume(
                                wellobj=wellobj, transfervolume=transfervolume
                            )
                            wellinfo.append(
                                [previousreactionobjs, wellobj, transfervolume]
                            )
                            transfervolume = 0.00
                        if wellvolumeavailable < transfervolume:
                            self.updateWellVolume(
                                wellobj=wellobj, transfervolume=wellvolumeavailable
                            )
                            wellinfo.append(
                                [previousreactionobjs, wellobj, wellvolumeavailable]
                            )
                            transfervolume = transfervolume - wellvolumeavailable
            except Exception as e:
                logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
                print(e)
                print(smiles, solvent, concentration)

        return wellinfo

    def checkVolumeClose(self, volume1: float, volume2: float) -> bool:
        """Checks if two volumes are almost the same value"""
        checkclose = math.isclose(volume1, volume2, rel_tol=0.001)
        if checkclose:
            return True
        else:
            return False

    def checkNextReactionsAddActions(
        self, reactionobj: Reaction, productsmiles: str
    ) -> list:
        """Checks if there are any reaction objects following the reaction in a method.
           If there is, checks if any of the proceeding reaction add actions match
           the reaction product's SMILES

        Parameters
        ----------
        reactionobj: Reaction
            The Django reaction model object to search for it's product SMILES
            matching any add actions needing the product as a reactant in the
            proceeding reactions
        productsmiles: str
            The SMILES of the reaction's product

        Returns
        -------
        addactionsmatches: list
            The Django add action model objects that require the reaction product
            as an input reactant
        """
        reactionqueryset = self.getReactionQuerySet(method_id=reactionobj.method_id.id)
        nextreactionqueryset = self.getNextObjEntries(
            queryset=reactionqueryset, obj=reactionobj
        )

        addactionsmatches = []
        for reactionobj in nextreactionqueryset:
            addactionmatch = self.getAddActionQuerySet(
                reaction_ids=[reactionobj.id],
                actionsessiontype="reaction",
            ).filter(smiles=productsmiles)
            if addactionmatch:
                addactionsmatches.append(addactionmatch[0])
        return addactionsmatches

    def getWellObj(self, reaction_id: int, welltype: str) -> Well:
        """Find the reaction plate well

        Parameters
        ----------
        reaction_id: int
            The reaction's id linked to the well on the reaction plate
        welltype: str
            The type of well eg.reaction, workup, lcms

        Returns
        -------
        wellobj: Well
            The well used in the reaction
        """
        productsmiles = self.getProductSmiles(reaction_id=reaction_id)
        # print(productsmiles, reaction_id, self.otsession_id, welltype)

        wellobj = Well.objects.get(
            otsession_id=self.otsession_id,
            reaction_id=reaction_id,
            type=welltype,
            smiles=productsmiles,
        )
        return wellobj

    def getWellVolumeAvailable(self, wellobj: Well) -> float:
        """Get volume of well available

        Parameters
        ----------
        wellobj: Well
            The well to check it's well volume availability
        """
        plateid = wellobj.plate_id.id
        maxwellvolume = self.getMaxWellVolume(plateid=plateid)
        deadvolume = self.getDeadVolume(maxwellvolume=maxwellvolume)
        wellvolume = wellobj.volume
        wellvolumeavailable = wellvolume - deadvolume
        self.updateWellAvailable(
            wellvolumeavailable=wellvolumeavailable, wellobj=wellobj
        )
        return wellvolumeavailable

    def updateWellAvailable(self, wellvolumeavailable: float, wellobj: Well):
        if wellvolumeavailable < 0:
            wellobj.available = False
            wellobj.save()

    def updateWellVolume(self, wellobj: Well, transfervolume: float):
        """Updates the volume available in a well after a transfer

        Parameters
        ---------
        wellobj: Well
            The well to update it's volume
        transfervolume: float
            The volume transferred from the well
        """
        wellobj.volume = wellobj.volume - transfervolume
        wellobj.save()

    def updateReactantsForNextStep(self, columnobj: Column):
        """Updates well objects related to a column
        to have reactant for next step set to True
        """
        columnwellqueryset = Well.objects.filter(column_id=columnobj)
        if columnwellqueryset:
            for columnwellobj in columnwellqueryset:
                clonewellqueryset = Well.objects.filter(id=columnwellobj.clonewellid)
                if clonewellqueryset:
                    clonewellobj = clonewellqueryset[0]
                    clonewellobj.reactantfornextstep = False
                    clonewellobj.save()
                columnwellobj.reactantfornextstep = True
                columnwellobj.save()

    def updateReactantsIsNotForNextStep(self, columnobj: Column):
        """Updates well objects related to a cloumn to have reactant for
        next step set to False
        """
        columnwellqueryset = Well.objects.filter(column_id=columnobj)
        if columnwellqueryset:
            for columnwellobj in columnwellqueryset:
                clonewellqueryset = Well.objects.filter(id=columnwellobj.clonewellid)
                if clonewellqueryset:
                    clonewellobj = clonewellqueryset[0]
                    clonewellobj.reactantfornextstep = False
                    clonewellobj.save()
                columnwellobj.reactantfornextstep = False
                columnwellobj.save()

    def updateReactantForNextStep(self, wellobj: Well):
        """Updates well object to have reactant for next step
        set to True
        """
        clonewellqueryset = Well.objects.filter(id=wellobj.clonewellid)
        if clonewellqueryset:
            clonewellobj = clonewellqueryset[0]
            clonewellobj.reactantfornextstep = False
            clonewellobj.save()
        wellobj.reactantfornextstep = True
        wellobj.save()

    def updateReactantIsNotForNextStep(self, wellobj: Well):
        """Updates well object to have reactant for next step
        set to False
        """
        clonewellqueryset = Well.objects.filter(id=wellobj.clonewellid)
        if clonewellqueryset:
            clonewellobj = clonewellqueryset[0]
            clonewellobj.reactantfornextstep = False
            clonewellobj.save()
        wellobj.reactantfornextstep = False
        wellobj.save()

    def getMaxWellVolume(self, plateid: int) -> float:
        """Get the maximum volume of a plate's wells

        Parameters
        ----------
        plateid: id
            The plate id to get the maximum well volume of

        Returns
        -------
        maxwellvolume: float
            The plate's maximum well volume
        """
        plateobj = self.getPlateObj(plateid)
        maxwellvolume = plateobj.maxwellvolume
        return maxwellvolume

    def getNextObjEntries(self, queryset: QuerySet, obj: DjangoObjectType) -> QuerySet:
        """Finds all proceeding Django model object relative to the Django model
           object in a queryset

        Parameters
        ----------
        queryset: QuerySet
            The queryset to search for proceeding entries
        obj: DjangoObjectType
            The object that you want to find all proceeding object entries relative to

        Returns
        -------
        nextqueryset: QuerySet
            The proceeding Django model objects as a queryset
        """
        nextqueryset = queryset.filter(pk__gt=obj.pk).order_by("pk")
        return nextqueryset

    def getDeadVolume(self, maxwellvolume: float) -> float:
        deadvolume = maxwellvolume * 0.05
        return deadvolume

    def getUniqueReactionClasses(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """Set of unique reaction classes

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to get unique list of reaction classes for

        Returns
        ------
        reactionclasses: list
            The set of reactionclasses
        """
        reactionclasses = (
            reactionqueryset.values_list("reactionclass", flat=True)
            .order_by("reactionclass")
            .distinct()
        )
        return reactionclasses

    def getGroupedReactionByClass(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """Group reactions by reaction class

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to to group by reaction class

        Returns
        -------
        groupedreactionbyclassquerysets: list
            The list of sublists of reaction querysets grouped by reaction class
        """
        reactionclasses = self.getUniqueReactionClasses(
            reactionqueryset=reactionqueryset
        )
        groupedreactionbyclassquerysets = []

        for reactionclass in reactionclasses:
            reactionbyclassqueryset = (
                reactionqueryset.filter(reactionclass=reactionclass)
                .distinct()
                .order_by("reactionclass")
            )
            if reactionbyclassqueryset:
                groupedreactionbyclassquerysets.append(reactionbyclassqueryset)

        return groupedreactionbyclassquerysets

    def getColumns(self, columntype: str, reactionclass: str) -> QuerySet[Column]:
        """Get columns related to a reactions of a partiular reaction class

        Parameters
        ----------
        columntype: str
            The type of column eg. reaction, workup1, workup2 etc
        reactionclass: str
            The reaction class to find columns for

        Returns
        -------
        columnqueryset: QuerySet[Columns]
            The columns containing wells where reactions of the same class
            were executed
        """
        columnqueryset = Column.objects.filter(
            otsession_id=self.otsession_id, type=columntype, reactionclass=reactionclass
        ).order_by("index")
        return columnqueryset

    def setupScript(self):
        """Writes header information for an OT script"""
        script = open(self.filepath, "w")
        script.write("from opentrons import protocol_api\n")
        script.write(
            "# "
            + str(self.protocolname)
            + str('" produced by XChem Car (https://car.xchem.diamond.ac.uk)')
        )
        script.write("\n# metadata")
        script.write(
            "\nmetadata = {'protocolName': '"
            + str(self.protocolname)
            + "','apiLevel': '"
            + str(self.apiLevel)
            + "'}\n"
        )
        script.write("\ndef run(protocol: protocol_api.ProtocolContext):\n")

        script.close()

    def setupPlates(self):
        """Writes the plate setup instructions for an OT script"""
        script = open(self.filepath, "a")
        script.write("\n\t# labware")
        for plateobj in self.platequeryset:
            platename = plateobj.name
            labware = plateobj.labware
            plateindex = plateobj.index
            script.write(
                f"\n\t{platename} = protocol.load_labware('{labware}', '{plateindex}')"
            )

        script.close()

    def setupTipRacks(self):
        """Writes the tipracks setup instructions for an OT script"""
        script = open(self.filepath, "a")
        for tiprackobj in self.tiprackqueryset:
            name = tiprackobj.name
            labware = tiprackobj.labware
            index = tiprackobj.index
            script.write(f"\n\t{name} = protocol.load_labware('{labware}', '{index}')")

        script.close()

    def setupPipettes(self):
        """Writes the pipette setup instructions for an OT script"""

        script = open(self.filepath, "a")
        script.write("\n\n\t# pipettes\n")
        script.write(
            "\t"
            + str(self.pipetteobj.name)
            + " = protocol.load_instrument('"
            + str(self.pipetteobj.labware)
            + "', '"
            + str(self.pipetteobj.position)
            + "', tip_racks=["
            + ",".join([tiprackobj.name for tiprackobj in self.tiprackqueryset])
            + "])\n"
        )

        script.close()

    def writeCommand(self, comandString: str):
        """Writes the command for an OT script"""
        script = open(self.filepath, "a")
        if type(comandString) == str:
            script.write("\t" + str(comandString) + "\n")
        elif type(comandString) == list:
            for command in comandString:
                script.write("\t" + str(command) + "\n")

        script.close()

    def mixWell(self, wellindex: int, nomixes: int, plate: str, volumetomix: float):
        """Prepares mixing commmand instruction for a mixing action

        Parameters
        ----------
        wellindex: int
            The index of the well to be mixed
        nomixes: int
            The number of mixes needed
        plate: str
            The plate linked to the well being mixed
        volumetomix: float
            The volume (ul) to be aspirated/dipsensed for the mxing action
        """
        humanread = f"Mixing contents of plate: {plate} at well index: {wellindex}"

        instruction = [
            "\n\t# " + str(humanread),
            self.pipettename
            + f".mix({nomixes}, {volumetomix}, {plate}.wells()[{wellindex}])",
        ]

        self.writeCommand(instruction)

    def mixColumn(self, columnindex: int, nomixes: int, plate: str, volumetomix: float):
        """Prepares mixing commmand instruction for a mixing action
        on a column in a plate

        Parameters
        ----------
        columnindex: int
            The index of the column to be mixed
        nomixes: int
            The number of mixes needed
        plate: str
            The plate linked to the column being mixed
        volumetomix: float
            The volume (ul) to be aspirated/dipsensed for the mxing action
        """
        humanread = f"Mixing contents of plate: {plate} at column index: {columnindex}"

        instruction = [
            "\n\t# " + str(humanread),
            self.pipettename
            + f".mix({nomixes}, {volumetomix}, {plate}[{columnindex}])",
        ]

        self.writeCommand(instruction)

    def pickUpTip(self):
        """Prepares a pick up tip action"""
        humanread = f"Picking up a new tip"

        instruction = [
            "\n\t# " + str(humanread),
            "if not {}.has_tip:\n".format(self.pipettename),
            "\t{}.pick_up_tip()".format(self.pipettename),
        ]

        self.writeCommand(instruction)

    def dropTip(self):
        """Prepares a drop a tip to waste"""
        humanread = f"Sending tip to waste"

        instruction = [
            "\n\t# " + str(humanread),
            "if {}.has_tip:\n".format(self.pipettename),
            "\t{}.drop_tip()".format(self.pipettename),
        ]

        self.writeCommand(instruction)

    def delayProtocol(self, delay: int):
        """Delays protocol from executing next operation

        Parameters
        ----------
        delay: int
            The delay in seconds

        """
        humanread = f"Delaying protocol operation"

        instruction = [
            "\n\t# " + str(humanread),
            "protocol.delay(seconds={})".format(delay),
        ]

        self.writeCommand(instruction)

    def transferFluidSingle(
        self,
        aspirateplatename: str,
        dispenseplatename: str,
        aspiratewellindex: int,
        dispensewellindex: int,
        transvolume: float,
        aspirateheight: int = 0.1,
        dispenseheight: int = -5,
        airgap: int = 15,
        transfertype: str = "standard",
    ):
        """Prepares the transfer commmand instruction for a tranfer action
           using a single well index

        Parameters
        ----------
        aspirateplatename: str
            The name of the plate the transfer is orignating from
        dispenseplatename: str
            The name of the destination plate
        aspiratewellindex: int
            The index of the well the transfer is orginating from
        dispensewellindex: int
            The index of the destination well
        transvolume: float
            The volume being transferred
        aspirateheight: int
            The height (mm) the pipette tip will aspirate from.
            Set to 2 mm.
        dispenseheight: int
            The hieght (mm) at which the pipette tip will dispense.
            Set to 5 mm above the bottom of the well.
        transfertype: str
            The tpye of transfer eg. reaction, workup
        airgap: int
            The airgap (uL) setting for the pipette to draw air in after aspirating
            to prevent solvent pearling
        newtip:
            Set to "never" to deal with pick up and drop tips built into protocol
        """

        humanread = f"transfertype - {transfertype} - transfer - {transvolume:.1f}ul from {aspiratewellindex} to {dispensewellindex}"
        instruction = [
            "\n\t# " + str(humanread),
            self.pipettename
            + f".transfer({transvolume}, {aspirateplatename}.wells()[{aspiratewellindex}].bottom({aspirateheight}), {dispenseplatename}.wells()[{dispensewellindex}].top({dispenseheight}), air_gap = {airgap}, new_tip='never', blow_out=True, blowout_location='destination well')",
        ]
        self.writeCommand(instruction)

    def transferFluidMulti(
        self,
        aspirateplatename: str,
        dispenseplatename: str,
        aspiratecolumnindex: int,
        dispensecolumnindex: int,
        transvolume: float,
        aspirateheight: int = 0.1,
        dispenseheight: int = -5,
        airgap: int = 15,
        transfertype: str = "standard",
    ):
        """Prepares the mutli-pipette transfer commmand instruction for a tranfer action
           using a column index

        Parameters
        ----------
        aspirateplatename: str
            The name of the plate the transfer is orignating from
        dispenseplatename: str
            The name of the destination plate
        aspiratecolumnindex: int
            The index of the column the transfer is orginating from
        dispensecolumnindex: int
            The index of the destination column
        transvolume: float
            The volume being transferred
        aspirateheight: int
            The height (mm) the pipette tip will aspirate from.
            Set to 2 mm.
        dispenseheight: int
            The hieght (mm) at which the pipette tip will dispense.
            Set to 5 mm above the bottom of the well.
        transfertype: str
            The tpye of transfer eg. reaction, workup
        airgap: int
            The airgap (uL) setting for the pipette to draw air in after aspirating
            to prevent solvent pearling
        newtip:
            Set to "never" to deal with pick up and drop tips built into protocol
        """
        print("Transferring multi")

        humanread = f"transfertype - {transfertype} - transfer - {transvolume:.1f}ul from {aspiratecolumnindex} column to {dispensecolumnindex} column"
        instruction = [
            "\n\t# " + str(humanread),
            self.pipettename
            + f".transfer({transvolume}, {aspirateplatename}[{aspiratecolumnindex}].bottom({aspirateheight}), {dispenseplatename}[{dispensecolumnindex}].top({dispenseheight}), air_gap = {airgap}, new_tip='never', blow_out=True, blowout_location='destination well')",
        ]
        self.writeCommand(instruction)

    def calculateAspirateHeight(self, bottomlayervolume: float, labware: str) -> float:
        """Calculate the height at wich to extract the top layer based
        on volume occupying the bottom layer

        Parameters
        ----------
        bottomlayervolume: float
            The volume occupying the bottom layer
        platetype: str
            The type of plate being used for the extraction

        Returns
        -------
        aspirateheight: float
            The height (mm) from the bottom of the plate,
            the pipette tip should aspirate from
        """
        aspirateheightconvesion_m = labware_plates[labware][
            "aspirateheightconversion-m"
        ]
        aspirateheightconvesion_c = labware_plates[labware][
            "aspirateheightconversion-c"
        ]
        bottomlayerheight = (
            aspirateheightconvesion_m * bottomlayervolume
        ) + aspirateheightconvesion_c
        aspirateheight = (bottomlayerheight * 0.15) + bottomlayerheight
        return aspirateheight

    def writeReactionActions(self, actionsessionqueryset: QuerySet[ActionSession]):
        for actionsessionobj in actionsessionqueryset:
            sessionnumber = actionsessionobj.sessionnumber
            reactionobj = Reaction.objects.get(id=actionsessionobj.reaction_id.id)
            reaction_id = reactionobj.id
            reactionclass = reactionobj.reactionclass
            recipetype = reactionobj.recipetype
            intramolecular = reactionobj.intramolecular
            if intramolecular:
                reactionactionsearch = "intramolecular"
            if not intramolecular:
                reactionactionsearch = "intermolecular"

            actionsessions = encoded_recipes[reactionclass]["recipes"][recipetype][
                "actionsessions"
            ]
            reactionactions = [
                actionsession[reactionactionsearch]["actions"]
                for actionsession in actionsessions
                if actionsession["type"] == "reaction"
                and actionsession["sessionnumber"] == sessionnumber
            ][0]

            for index, reactionaction in enumerate(reactionactions):
                actiontype = reactionaction["type"]
                actionnumber = reactionaction["actionnumber"]
                if actiontype == "add":
                    toplatetype = reactionaction["content"]["plates"]["toplatetype"]
                    addactionobj = AddAction.objects.get(
                        actionsession_id=actionsessionobj,
                        reaction_id=reactionobj,
                        number=actionnumber,
                    )
                    smiles = addactionobj.smiles
                    solvent = addactionobj.solvent
                    transfervolume = addactionobj.volume
                    concentration = addactionobj.concentration

                    fromwellinfo = self.findStartingPlateWellObj(
                        reaction_id=reaction_id,
                        smiles=smiles,
                        solvent=solvent,
                        concentration=concentration,
                        transfervolume=transfervolume,
                    )

                    for wellinfo in fromwellinfo:
                        previousreactionobjs = wellinfo[0]
                        fromwellobj = wellinfo[1]
                        transfervolume = wellinfo[2]

                        if previousreactionobjs:
                            fromsolventwellinfo = self.findSolventPlateWellObj(
                                solvent=solvent,
                                transfervolume=transfervolume,
                            )

                            self.pickUpTip()
                            for solventwellinfo in fromsolventwellinfo:
                                fromsolventwellobj = solventwellinfo[0]
                                transfervolume = solventwellinfo[1]
                                towellobj = fromwellobj
                                fromplateobj = self.getPlateObj(
                                    plateid=fromsolventwellobj.plate_id.id
                                )
                                toplateobj = self.getPlateObj(
                                    plateid=towellobj.plate_id.id
                                )

                                aspirateplatename = fromplateobj.name
                                dispenseplatename = toplateobj.name
                                aspiratewellindex = fromsolventwellobj.index
                                dispensewellindex = towellobj.index

                                self.transferFluidSingle(
                                    aspirateplatename=aspirateplatename,
                                    dispenseplatename=dispenseplatename,
                                    aspiratewellindex=aspiratewellindex,
                                    dispensewellindex=dispensewellindex,
                                    transvolume=transfervolume,
                                    transfertype="dilution",
                                )

                        towellobj = self.getWellObj(
                            reaction_id=reaction_id,
                            welltype=toplatetype,
                        )
                        fromplateobj = self.getPlateObj(plateid=fromwellobj.plate_id.id)
                        toplateobj = self.getPlateObj(plateid=towellobj.plate_id.id)

                        aspirateplatename = fromplateobj.name
                        dispenseplatename = toplateobj.name
                        aspiratewellindex = fromwellobj.index
                        dispensewellindex = towellobj.index

                        self.pickUpTip()
                        self.transferFluidSingle(
                            aspirateplatename=aspirateplatename,
                            dispenseplatename=dispenseplatename,
                            aspiratewellindex=aspiratewellindex,
                            dispensewellindex=dispensewellindex,
                            transvolume=transfervolume,
                        )
                        if index + 1 == len(reactionactions):
                            self.dropTip()
                        if index + 1 < len(reactionactions):
                            if reactionactions[index + 1]["type"] != "mix":
                                self.dropTip()
                        self.updateReactantForNextStep(wellobj=towellobj)
                        self.updateReactantIsNotForNextStep(wellobj=fromwellobj)

                if actiontype == "mix":
                    mixactionobj = MixAction.objects.get(
                        actionsession_id=actionsessionobj,
                        reaction_id=reactionobj,
                        number=actionnumber,
                    )
                    platetype = mixactionobj.platetype
                    repetitions = mixactionobj.repetitions

                    mixwellobj = self.getWellObj(
                        reaction_id=reaction_id, welltype=platetype
                    )
                    mixplateobj = self.getPlateObj(plateid=mixwellobj.plate_id.id)
                    mixwellindex = mixwellobj.index
                    mixplatename = mixplateobj.name

                    self.mixWell(
                        wellindex=mixwellindex,
                        nomixes=repetitions,
                        plate=mixplatename,
                        volumetomix=transfervolume,
                    )
                    self.dropTip()

    def writeWorkUpActions(self, actionsessionqueryset: QuerySet[ActionSession]):
        sessionnumber = actionsessionqueryset.values_list(
            "sessionnumber", flat=True
        ).distinct()[0]
        reaction_ids = actionsessionqueryset.values_list(
            "reaction_id", flat=True
        ).order_by("reaction_id")
        reactionqueryset = self.getReactionQuerySet(reaction_ids=reaction_ids)
        groupedreactionclassquerysets = self.getGroupedReactionByClass(
            reactionqueryset=reactionqueryset
        )
        for groupreactionclassqueryset in groupedreactionclassquerysets:
            actionsessiontype = "workup"
            reactionclass = groupreactionclassqueryset.values_list(
                "reactionclass", flat=True
            ).distinct()[0]
            recipetype = "standard"  # Introduce loop here to group by recipe type
            actionsessions = encoded_recipes[reactionclass]["recipes"][recipetype][
                "actionsessions"
            ]
            workupactions = [
                actionsession["actions"]
                for actionsession in actionsessions
                if actionsession["type"] == "workup"
                and actionsession["sessionnumber"] == sessionnumber
            ][0]
            for index, workupaction in enumerate(workupactions):
                actiontype = workupaction["type"]
                actionnumber = workupaction["actionnumber"]
                if actiontype == "add":
                    addactionqueryset = self.getAddActionQuerySet(
                        reaction_ids=groupreactionclassqueryset,
                        actionsessiontype=actionsessiontype,
                        actionnumber=actionnumber,
                    )
                    transfervolume = addactionqueryset.values_list(
                        "volume", flat=True
                    ).distinct()[0]
                    solvent = addactionqueryset.values_list(
                        "solvent", flat=True
                    ).distinct()[0]
                    toplatetype = workupaction["content"]["plates"]["toplatetype"]
                    tocolumnqueryset = self.getColumnQuerySet(
                        columntype=toplatetype, reactionclass=reactionclass
                    )
                    self.pickUpTip()
                    for topcolumnobj in tocolumnqueryset:
                        toplateobj = topcolumnobj.plate_id
                        dispenseplatename = toplateobj.name
                        dispensecolumnindex = topcolumnobj.index

                        fromsolventwellinfo = self.findSolventPlateWellObj(
                            solvent=solvent,
                            transfervolume=transfervolume,
                        )
                        for solventwellinfo in fromsolventwellinfo:
                            fromsolventwellobj = solventwellinfo[0]
                            transfervolume = solventwellinfo[1]

                            fromplateobj = self.getPlateObj(
                                plateid=fromsolventwellobj.plate_id.id
                            )
                            aspirateplatename = fromplateobj.name
                            aspiratecolumnindex = fromsolventwellobj.index

                            self.transferFluidMulti(
                                aspirateplatename=aspirateplatename,
                                dispenseplatename=dispenseplatename,
                                aspiratecolumnindex=aspiratecolumnindex,
                                dispensecolumnindex=dispensecolumnindex,
                                transvolume=transfervolume,
                                transfertype="workup",
                            )
                    self.dropTip()

                    if index + 1 == len(workupactions):
                        self.dropTip()
                    if index + 1 < len(workupactions):
                        if workupactions[index + 1]["type"] != "mix":
                            self.dropTip()

                if actiontype == "extract":
                    fromplatetype = workupaction["content"]["plates"]["fromplatetype"]
                    toplatetype = workupaction["content"]["plates"]["toplatetype"]
                    extractactionqueyset = self.getExtractActionQuerySet(
                        reaction_ids=groupreactionclassqueryset,
                        actionsessiontype=actionsessiontype,
                        actionnumber=actionnumber,
                    )
                    transfervolume = extractactionqueyset.values_list(
                        "volume", flat=True
                    ).distinct()[0]
                    solvent = extractactionqueyset.values_list(
                        "solvent", flat=True
                    ).distinct()[0]
                    bottomlayervolume = extractactionqueyset.values_list(
                        "bottomlayervolume", flat=True
                    ).distinct()[0]
                    fromcolumnqueryset = self.getColumnQuerySet(
                        columntype=fromplatetype, reactionclass=reactionclass
                    )
                    tocolumnqueryset = self.getColumnQuerySet(
                        columntype=toplatetype, reactionclass=reactionclass
                    )
                    aspirateheight = None
                    for fromcolumnobj, tocolumnobj in zip(
                        fromcolumnqueryset, tocolumnqueryset
                    ):
                        fromplateobj = fromcolumnobj.plate_id
                        aspirateplatename = fromplateobj.name
                        aspiratecolumnindex = fromcolumnobj.index
                        toplateobj = tocolumnobj.plate_id
                        dispenseplatename = toplateobj.name
                        dispensecolumnindex = tocolumnobj.index
                        if not aspirateheight:
                            if bottomlayervolume:
                                aspirateheight = self.calculateAspirateHeight(
                                    labware=fromplateobj.labware,
                                    bottomlayervolume=bottomlayervolume,
                                )
                            else:
                                aspirateheight = 0.1
                        self.delayProtocol(delay=5)
                        self.pickUpTip()
                        self.transferFluidMulti(
                            aspirateplatename=aspirateplatename,
                            dispenseplatename=dispenseplatename,
                            aspiratecolumnindex=aspiratecolumnindex,
                            dispensecolumnindex=dispensecolumnindex,
                            transvolume=transfervolume,
                            transfertype="workup",
                            aspirateheight=aspirateheight,
                        )
                        self.dropTip()
                        self.updateReactantsForNextStep(columnobj=tocolumnobj)
                        self.updateReactantsIsNotForNextStep(columnobj=fromcolumnobj)

                    if index + 1 == len(workupactions):
                        self.dropTip()
                    if index + 1 < len(workupactions):
                        if workupactions[index + 1]["type"] != "mix":
                            self.dropTip()

                if actiontype == "mix":
                    mixactionqueyset = self.getMixActionQuerySet(
                        reaction_ids=groupreactionclassqueryset,
                        actionsessiontype=actionsessiontype,
                        actionnumber=actionnumber,
                    )
                    platetype = mixactionqueyset.values_list(
                        "platetype", flat=True
                    ).distinct()[0]
                    repetitions = mixactionqueyset.values_list(
                        "repetitions", flat=True
                    ).distinct()[0]
                    mixvolume = mixactionqueyset.values_list(
                        "volume", flat=True
                    ).distinct()[0]

                    mixcolumnqueryset = self.getColumnQuerySet(
                        columntype=platetype, reactionclass=reactionclass
                    )

                    for mixcolumnobj in mixcolumnqueryset:
                        self.pickUpTip()
                        mixplateobj = mixcolumnobj.plate_id
                        mixplatename = mixplateobj.name
                        mixcolumnindex = mixcolumnobj.index
                        self.mixColumn(
                            columnindex=mixcolumnindex,
                            nomixes=repetitions,
                            plate=mixplatename,
                            volumetomix=mixvolume,
                        )
                        self.dropTip()

    def writeAnalyseActions(self, actionsessionqueryset: QuerySet[ActionSession]):
        sessionnumber = actionsessionqueryset.values_list(
            "sessionnumber", flat=True
        ).distinct()[0]
        reaction_ids = actionsessionqueryset.values_list(
            "reaction_id", flat=True
        ).order_by("reaction_id")
        reactionqueryset = self.getReactionQuerySet(reaction_ids=reaction_ids)
        groupedreactionclassquerysets = self.getGroupedReactionByClass(
            reactionqueryset=reactionqueryset
        )

        for groupreactionclassqueryset in groupedreactionclassquerysets:
            actionsessiontype = "analyse"
            reactionclass = groupreactionclassqueryset.values_list(
                "reactionclass", flat=True
            ).distinct()[0]
            recipetype = "standard"  # Introduce loop here to group by recipe type
            actionsessions = encoded_recipes[reactionclass]["recipes"][recipetype][
                "actionsessions"
            ]
            analyseactions = [
                actionsession["actions"]
                for actionsession in actionsessions
                if actionsession["type"] == "analyse"
                and actionsession["sessionnumber"] == sessionnumber
            ][0]

            for index, analyseaction in enumerate(analyseactions):
                actiontype = analyseaction["type"]
                actionnumber = analyseaction["actionnumber"]
                if actiontype == "add":
                    fromplatetype = analyseaction["content"]["plates"]["fromplatetype"]
                    toplatetype = analyseaction["content"]["plates"]["toplatetype"]
                    addactionqueryset = self.getAddActionQuerySet(
                        reaction_ids=groupreactionclassqueryset,
                        actionsessiontype=actionsessiontype,
                        actionnumber=actionnumber,
                    )
                    transfervolume = addactionqueryset.values_list(
                        "volume", flat=True
                    ).distinct()[0]
                    solvent = addactionqueryset.values_list(
                        "solvent", flat=True
                    ).distinct()[0]
                    tocolumnqueryset = self.getColumnQuerySet(
                        columntype=toplatetype, reactionclass=reactionclass
                    )

                    if fromplatetype == "solvent":
                        print("Adding solvent")
                        self.pickUpTip()
                        print("The to column queryset is: {}".format(tocolumnqueryset))
                        for topcolumnobj in tocolumnqueryset:
                            toplateobj = topcolumnobj.plate_id
                            dispenseplatename = toplateobj.name
                            dispensecolumnindex = topcolumnobj.index

                            fromsolventwellinfo = self.findSolventPlateWellObj(
                                solvent=solvent,
                                transfervolume=transfervolume,
                            )
                            print(
                                "The from solvent well info is: {}".format(
                                    fromsolventwellinfo
                                )
                            )
                            for solventwellinfo in fromsolventwellinfo:
                                fromsolventwellobj = solventwellinfo[0]
                                transfervolume = solventwellinfo[1]

                                fromplateobj = self.getPlateObj(
                                    plateid=fromsolventwellobj.plate_id.id
                                )
                                aspirateplatename = fromplateobj.name
                                aspiratecolumnindex = fromsolventwellobj.index

                                self.transferFluidMulti(
                                    aspirateplatename=aspirateplatename,
                                    dispenseplatename=dispenseplatename,
                                    aspiratecolumnindex=aspiratecolumnindex,
                                    dispensecolumnindex=dispensecolumnindex,
                                    transvolume=transfervolume,
                                    transfertype="analyse",
                                )
                        self.dropTip()

                    if fromplatetype in ["reaction", "workup1", "workup2", "workup3"]:
                        fromcolumnqueryset = self.getColumnQuerySet(
                            columntype=fromplatetype, reactionclass=reactionclass
                        )
                        for fromcolumnobj, tocolumnobj in zip(
                            fromcolumnqueryset, tocolumnqueryset
                        ):
                            fromplateobj = fromcolumnobj.plate_id
                            aspirateplatename = fromplateobj.name
                            aspiratecolumnindex = fromcolumnobj.index
                            toplateobj = tocolumnobj.plate_id
                            dispenseplatename = toplateobj.name
                            dispensecolumnindex = tocolumnobj.index
                            self.delayProtocol(delay=5)
                            self.pickUpTip()
                            self.transferFluidMulti(
                                aspirateplatename=aspirateplatename,
                                dispenseplatename=dispenseplatename,
                                aspiratecolumnindex=aspiratecolumnindex,
                                dispensecolumnindex=dispensecolumnindex,
                                transvolume=transfervolume,
                                transfertype="analyse",
                            )
                            self.dropTip()

                    if index + 1 == len(analyseactions):
                        self.dropTip()
                    if index + 1 < len(analyseactions):
                        if analyseactions[index + 1]["type"] != "mix":
                            self.dropTip()

                if actiontype == "mix":
                    mixactionqueyset = self.getMixActionQuerySet(
                        reaction_ids=groupreactionclassqueryset,
                        actionsessiontype=actionsessiontype,
                        actionnumber=actionnumber,
                    )
                    platetype = mixactionqueyset.values_list(
                        "platetype", flat=True
                    ).distinct()[0]
                    repetitions = mixactionqueyset.values_list(
                        "repetitions", flat=True
                    ).distinct()[0]
                    mixvolume = mixactionqueyset.values_list(
                        "volume", flat=True
                    ).distinct()[0]

                    mixcolumnqueryset = self.getColumnQuerySet(
                        columntype=platetype, reactionclass=reactionclass
                    )

                    for mixcolumnobj in mixcolumnqueryset:
                        self.pickUpTip()
                        mixplateobj = mixcolumnobj.plate_id
                        mixplatename = mixplateobj.name
                        mixcolumnindex = mixcolumnobj.index
                        self.mixColumn(
                            columnindex=mixcolumnindex,
                            nomixes=repetitions,
                            plate=mixplatename,
                            volumetomix=mixvolume,
                        )
                        self.dropTip()
