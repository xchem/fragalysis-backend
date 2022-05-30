# this is vunrable to python injection by the lack of checking of metadata inputs
# this opens and closes files frequently, could be improved by creating string to hold the file data before writing to file once

#   NOTE: in this file "humanread" referes to the comments above each line/set of lines of ot code, human readable is a list of all the comments in format [oporator (human/ot), comment]
"""Create otWrite object"""
from __future__ import annotations
from django.core.files.storage import default_storage
from django.conf import settings
from django.db.models import QuerySet

import os

from graphene_django import DjangoObjectType

from car.models import AddAction, AnalyseAction, OTSession, Reaction

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

logger = logging.getLogger(__name__)


class otWrite(object):
    """ "
    Creates a otWrite object for generating an OT protocol
    script
    """

    def __init__(
        self,
        protocolname: str,
        otsessionobj: OTSession,
        apiLevel: str = "2.9",
        alladdactionsquerysetflat: list[AddAction] = None,
        allanalyseactionsqueryset: QuerySet[AnalyseAction] = None,
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
        alladdactionsquerysetflat: list[AddAction]
            The optional add actions that the protocol is being created for
        allanalyseactionsqueryset: QuerySet[AnalyseAction]
            The optional analyse actions the protocol is being created for
        """
        self.reactionstep = otsessionobj.reactionstep
        self.otsessionobj = otsessionobj
        self.otsessionid = otsessionobj.id
        self.otsessiontype = otsessionobj.sessiontype
        self.protocolname = protocolname
        self.apiLevel = apiLevel
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
            self.writeReactionSession(
                alladdactionsquerysetflat=alladdactionsquerysetflat,
            )
        if self.otsessiontype == "analyse":
            self.writeAnalyseSession(
                allanalyseactionsqueryset=allanalyseactionsqueryset
            )

    def writeAnalyseSession(self, allanalyseactionsqueryset: QuerySet[AnalyseAction]):
        self.allanalyseactionsqueryset = allanalyseactionsqueryset
        self.writeAnalyseActions()
        self.createOTScriptModel()

    def writeReactionSession(self, alladdactionsquerysetflat: list[AddAction]):
        self.alladdactionsquerysetflat = alladdactionsquerysetflat
        self.writeAddActions()
        self.createOTScriptModel()

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
        platequeryset = Plate.objects.filter(otsession_id=self.otsessionid).order_by(
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
            otsession_id=self.otsessionid
        ).order_by("id")
        return tipracksqueryset

    def getPipette(self) -> Pipette:
        """Get the pipette for an OT session

        Returns
        -------
        pipetteobj: Pipette
            The pipette used for an OT session
        """

        pipetteobj = Pipette.objects.get(otsession_id=self.otsessionid)
        return pipetteobj

    def getProductSmiles(self, reactionid: int) -> str:
        """Get the product smiles

        Parameters
        ----------
        reactionid: int
            The reaction id linked to a product

        Returns
        -------
        productobj.smiles: str
            The SMILES of the product of  reaction
        """
        productobj = Product.objects.filter(reaction_id=reactionid)[0]
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

    def getReactionQuerySet(self, method_id: int) -> QuerySet[Reaction]:
        """Get a  synthesis methods reactions

        Parameters
        ----------
        method_id: int
            The synthesis method's id to get reactions for

        Returns
        -------
        reactionqueryset: QuerySet[Reaction]
            The reactions of a synthesis method
        """
        reactionqueryset = Reaction.objects.filter(method_id=method_id)
        return reactionqueryset

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
                self.otsessionid,
            )
        )
        path = "tmp/" + filename
        filepath = str(os.path.join(settings.MEDIA_ROOT, path))
        return filepath, filename

    def createOTScriptModel(self):
        """Creates an OT Script model object"""
        otscriptobj = OTScript()
        otscriptobj.otsession_id = self.otsessionobj
        otscriptfile = open(self.filepath, "rb")
        otscriptfn = default_storage.save(
            "otscripts/{}.py".format(self.filename.strip(".txt")), otscriptfile
        )
        otscriptobj.otscript = otscriptfn
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
                otsession_id=self.otsessionid,
                solvent=solvent,
                available=True,
                type="dilution",
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
        reactionid: int,
        smiles: str,
        solvent: str,
        concentration: float,
        transfervolume: float,
    ) -> list:
        """Finds starting plate well for executing an add action

        Parameters
        ----------
        reactionid: int
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
            reaction_id=reactionid, smiles=smiles
        )
        wellinfo = []
        if previousreactionobjs:
            wellobj = Well.objects.get(
                otsession_id=self.otsessionid,
                reaction_id=previousreactionobjs[0],
                smiles=smiles,
                type="reaction",
            )
            wellinfo.append([previousreactionobjs, wellobj, transfervolume])
        else:
            try:
                wellobjects = Well.objects.filter(
                    otsession_id=self.otsessionid,
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

    def findReactionPlateWellObj(self, reactionid: int) -> Well:
        """Find the reaction plate well

        Parameters
        ----------
        reactionid: int
            The reaction's id linked to the well on the reaction plate

        Returns
        -------
        wellobj: Well
            The well used in the reaction
        """
        productsmiles = self.getProductSmiles(reactionid=reactionid)
        wellobj = Well.objects.get(
            otsession_id=self.otsessionid, reaction_id=reactionid, smiles=productsmiles
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

    def getDeadVolume(self, maxwellvolume: float) -> float:
        deadvolume = maxwellvolume * 0.05
        return deadvolume

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
            + f".mix({nomixes}, {volumetomix}, {plate}.wells([{wellindex}]))",
        ]

        self.writeCommand(instruction)

    def transferFluid(
        self,
        aspirateplatename: str,
        dispenseplatename: str,
        aspiratewellindex: int,
        dispensewellindex: int,
        transvolume: float,
        aspirateheight: int = 2,
        dispenseheight: int = -5,
        transfertype: str = "standard",
    ):
        """Prepares the transfer commmand instruction for a tranfer action

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
        aspirateheight: int = 2
            The height the pipette tip will aspirate from
        dispenseheight: int = -5
            The hieght at which the pipette tip will dispense
        transfertype: str = "standard",
        """

        humanread = f"transfertype - {transfertype} - transfer - {transvolume:.1f}ul from {aspiratewellindex} to {dispensewellindex}"

        instruction = [
            "\n\t# " + str(humanread),
            self.pipettename
            + f".transfer({transvolume}, {aspirateplatename}.wells()[{aspiratewellindex}].bottom({aspirateheight}), {dispenseplatename}.wells()[{dispensewellindex}].top({dispenseheight}), air_gap = 15)",
        ]

        self.writeCommand(instruction)

    def writeAddActions(self):
        for addaction in self.alladdactionsquerysetflat:
            transfervolume = addaction.volume
            solvent = addaction.solvent

            fromwellinfo = self.findStartingPlateWellObj(
                reactionid=addaction.reaction_id.id,
                smiles=addaction.smiles,
                solvent=solvent,
                concentration=addaction.concentration,
                transfervolume=transfervolume,
            )

            for wellinfo in fromwellinfo:
                previousreactionobjs = wellinfo[0]
                fromwellobj = wellinfo[1]
                transfervolume = wellinfo[2]

                if previousreactionobjs:
                    # Add dilution
                    fromsolventwellinfo = self.findSolventPlateWellObj(
                        solvent=solvent,
                        transfervolume=transfervolume,
                    )
                    for solventwellinfo in fromsolventwellinfo:
                        fromsolventwellobj = solventwellinfo[0]
                        transfervolume = solventwellinfo[1]
                        towellobj = fromwellobj
                        fromplateobj = self.getPlateObj(
                            plateid=fromsolventwellobj.plate_id.id
                        )
                        toplateobj = self.getPlateObj(plateid=towellobj.plate_id.id)

                        aspirateplatename = fromplateobj.name
                        dispenseplatename = toplateobj.name
                        aspiratewellindex = fromsolventwellobj.index
                        dispensewellindex = towellobj.index

                        self.transferFluid(
                            aspirateplatename=aspirateplatename,
                            dispenseplatename=dispenseplatename,
                            aspiratewellindex=aspiratewellindex,
                            dispensewellindex=dispensewellindex,
                            transvolume=transfervolume,
                            transfertype="dilution",
                        )

                    self.mixWell(
                        wellindex=dispensewellindex,
                        nomixes=3,
                        plate=dispenseplatename,
                        volumetomix=transfervolume,
                    )

                towellobj = self.findReactionPlateWellObj(
                    reactionid=addaction.reaction_id.id
                )
                fromplateobj = self.getPlateObj(plateid=fromwellobj.plate_id.id)
                toplateobj = self.getPlateObj(plateid=towellobj.plate_id.id)

                aspirateplatename = fromplateobj.name
                dispenseplatename = toplateobj.name
                aspiratewellindex = fromwellobj.index
                dispensewellindex = towellobj.index

                self.transferFluid(
                    aspirateplatename=aspirateplatename,
                    dispenseplatename=dispenseplatename,
                    aspiratewellindex=aspiratewellindex,
                    dispensewellindex=dispensewellindex,
                    transvolume=transfervolume,
                )

    def writeAnalyseActions(self):
        for analyseaction in self.allanalyseactionsqueryset:
            samplevolume = analyseaction.samplevolume
            analysesolvent = analyseaction.solvent
            solventvolume = analyseaction.solventvolume
            reaction_id = analyseaction.reaction_id.id
            productsmiles = self.getProductSmiles(reactionid=reaction_id)

            fromwellobj = Well.objects.get(
                otsession_id=self.otsessionid,
                reaction_id=reaction_id,
                type="reaction",
                smiles=productsmiles,
            )
            fromplateobj = self.getPlateObj(plateid=fromwellobj.plate_id.id)
            towellobj = Well.objects.get(
                otsession_id=self.otsessionid,
                reaction_id=reaction_id,
                type="analyse",
                smiles=productsmiles,
            )
            toplateobj = self.getPlateObj(plateid=towellobj.plate_id.id)

            aspirateplatename = fromplateobj.name
            dispenseplatename = toplateobj.name
            aspiratewellindex = fromwellobj.index
            dispensewellindex = towellobj.index

            self.transferFluid(
                aspirateplatename=aspirateplatename,
                dispenseplatename=dispenseplatename,
                aspiratewellindex=aspiratewellindex,
                dispensewellindex=dispensewellindex,
                transvolume=samplevolume,
            )

            fromsolventwellinfo = self.findSolventPlateWellObj(
                solvent=analysesolvent,
                transfervolume=solventvolume,
            )

            for solventwellinfo in fromsolventwellinfo:
                fromsolventwellobj = solventwellinfo[0]
                transfervolume = solventwellinfo[1]

                fromplateobj = self.getPlateObj(plateid=fromsolventwellobj.plate_id.id)
                aspirateplatename = fromplateobj.name
                aspiratewellindex = fromsolventwellobj.index

                self.transferFluid(
                    aspirateplatename=aspirateplatename,
                    dispenseplatename=dispenseplatename,
                    aspiratewellindex=aspiratewellindex,
                    dispensewellindex=dispensewellindex,
                    transvolume=transfervolume,
                    transfertype="dilution",
                )

            self.mixWell(
                wellindex=dispensewellindex,
                nomixes=3,
                plate=dispenseplatename,
                volumetomix=solventvolume,
            )
