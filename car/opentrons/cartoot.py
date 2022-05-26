"""Create OT session"""
from __future__ import annotations
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from statistics import median
from graphene_django import DjangoObjectType

from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
from pandas.core.frame import DataFrame


from car.models import (
    AnalyseAction,
    Batch,
    Reaction,
    Product,
    AddAction,
    OTSession,
    Deck,
    Plate,
    SolventPrep,
    Well,
    Pipette,
    TipRack,
    CompoundOrder,
)

import math
from .labwareavailable import labware_plates


class CreateOTSession(object):
    """
    Creates a StartOTSession object for generating a protocol
    from actions for reactionqueryset for executing reactions
    """

    def __init__(
        self,
        reactionstep: int,
        otbatchprotocolobj: DjangoObjectType,
        sessiontype: str,
        reactiongroupqueryset: list = None,
    ):
        self.reactionstep = reactionstep
        self.otbatchprotocolobj = otbatchprotocolobj
        self.otsessionqueryset = self.otbatchprotocolobj.otsessions.all()
        self.batchobj = Batch.objects.get(id=otbatchprotocolobj.batch_id_id)
        self.sessiontype = sessiontype
        if sessiontype == "reaction":
            self.createReactionSession(reactiongroupqueryset=reactiongroupqueryset)
        if sessiontype == "analyse":
            self.createAnalyseSession(reactiongroupqueryset=reactiongroupqueryset)

    def createAnalyseSession(self, reactiongroupqueryset):
        self.otsessionobj = self.createOTSessionModel()
        self.reactiongroupqueryset = reactiongroupqueryset
        self.reactionids = [
            reactionobj.id for reactionobj in self.reactiongroupqueryset
        ]
        self.allanalyseactionqueryset = self.getAnalyseActionQuerySet(
            reaction_ids=self.reactionids
        )
        self.productqueryset = self.getProductQuerySet(reaction_ids=self.reactionids)
        self.productsmiles = [productobj.smiles for productobj in self.productqueryset]

        self.analyseactionsdf = self.getAnalyseActionsDataFrame()
        self.roundedvolumes = self.getRoundedAnalyseActionVolumes(
            analyseactionqueryset=self.allanalyseactionqueryset
        )
        self.groupedanalysemethodobjs = self.getGroupedAnalyseMethods()
        self.deckobj = self.createDeckModel()
        self.numbertips = self.getNumberTips(queryset=self.allanalyseactionqueryset)
        self.tipracktype = self.getTipRackType(roundedvolumes=self.roundedvolumes)
        self.createTipRacks(tipracktype=self.tipracktype)
        self.pipettetype = self.getPipetteType(roundedvolumes=self.roundedvolumes)
        self.createPipetteModel()
        self.createSolventPlate(materialsdf=self.analyseactionsdf)
        self.inputplatequeryset = self.getPreviousOTSessionReactionPlates()
        if self.inputplatequeryset:
            self.cloneInputPlate()
        for analysegroup in self.groupedanalysemethodobjs:
            # May need to rethink and do seprate sessions (Deck space) or optimise tip usage
            wellsneeded = len(analysegroup)
            method = analysegroup[0].method
            self.platetype = self.getAnalysePlateType(
                wellsneeded=wellsneeded, method=method
            )
            self.createAnalysePlate(
                method=method,
                platetype=self.platetype,
                analyseactionqueryset=analysegroup,
            )

    def createReactionSession(self, reactiongroupqueryset):
        self.otsessionobj = self.createOTSessionModel()
        self.reactiongroupqueryset = reactiongroupqueryset
        self.groupedtemperaturereactionobjs = self.getGroupedTemperatureReactions()
        self.alladdactionqueryset = [
            self.getAddActions(reaction_id=reactionobj.id)
            for reactionobj in self.reactiongroupqueryset
        ]
        self.alladdactionquerysetflat = [
            item for sublist in self.alladdactionqueryset for item in sublist
        ]
        self.roundedvolumes = self.getRoundedAddActionVolumes(
            addactionqueryset=self.alladdactionquerysetflat
        )
        self.deckobj = self.createDeckModel()
        self.numbertips = self.getNumberTips(queryset=self.alladdactionquerysetflat)
        self.tipracktype = self.getTipRackType(roundedvolumes=self.roundedvolumes)
        self.createTipRacks(tipracktype=self.tipracktype)
        self.pipettetype = self.getPipetteType(roundedvolumes=self.roundedvolumes)
        self.addactionsdf = self.getAddActionsDataFrame()
        self.inputplatequeryset = self.getInputPlatesNeeded()
        if self.inputplatequeryset:
            self.cloneInputPlate()
        self.createPipetteModel()
        self.createReactionStartingPlate()
        self.solventmaterialsdf = self.getAddActionsMaterialDataFrame(
            productexists=True
        )
        self.createSolventPlate(materialsdf=self.solventmaterialsdf)
        self.createReactionPlate()
        self.startingreactionplatequeryset = self.getStartingReactionPlateQuerySet()

    def getAllPreviousOTSessionReactionPlates(self):
        """Get all input reaction plates for all previous otsessions
        of reaction type in a otbatchprotocol
        """
        if self.otsessionqueryset:
            otsessionsids = [
                otsession.id
                for otsession in self.otsessionqueryset
                if otsession.sessiontype == "reaction"
            ]
            allinputplatesqueryset = Plate.objects.filter(
                otsession_id__in=otsessionsids, type="reaction"
            )
            return allinputplatesqueryset
        else:
            return False

    def getPlateQuerySet(self, otsession_id: int):
        platequeryset = Plate.objects.filter(otsession_id=otsession_id)
        return platequeryset

    def getPreviousOTSessionReactionPlates(self):
        """Checks if previous Reaction Plates exist and if they do
        check if the Wells match the products of the input
        reactiongroupqueryset for the session
        """
        previousotsessionobj = self.getPreviousObjEntry(
            queryset=self.otsessionqueryset, obj=self.otsessionobj
        )
        inputplatesneeded = []
        if previousotsessionobj.sessiontype == "reaction":
            previousotsessionplates = self.getPlateQuerySet(
                otsession_id=previousotsessionobj.id
            )
            previousotsessioreactionplates = previousotsessionplates.filter(
                type="reaction"
            )
            for previousotsessionreactionplate in previousotsessioreactionplates:
                wellmatchqueryset = (
                    previousotsessionreactionplate.well_set.all()
                    .filter(
                        reaction_id__in=self.reactionids,
                        smiles__in=self.productsmiles,
                        type="reaction",
                    )
                    .distinct()
                )
                if wellmatchqueryset:
                    inputplatesneeded.append(previousotsessionreactionplate)

        return inputplatesneeded

    def getInputPlatesNeeded(self):
        """Checks if a plate previosuly created is needed for the current otsession
        reactions
        """
        inputplatesneeded = []
        allinputplatequerset = self.getAllPreviousOTSessionReactionPlates()
        methodids = [
            reactionobj.method_id for reactionobj in self.reactiongroupqueryset
        ]
        allreactantsmiles = [
            addaction.smiles for addaction in self.alladdactionquerysetflat
        ]
        if allinputplatequerset:
            for inputplateobj in allinputplatequerset:
                wellmatchqueryset = (
                    inputplateobj.well_set.all()
                    .filter(
                        method_id__in=methodids,
                        reactantfornextstep=True,
                        smiles__in=allreactantsmiles,
                        type="reaction",
                    )
                    .distinct()
                )
                if wellmatchqueryset:
                    inputplatesneeded.append(inputplateobj)
        return inputplatesneeded

    def getPreviousObjEntry(self, queryset: list, obj: DjangoObjectType):
        """Finds all previous objects relative to obj of queryset"""
        previousobj = queryset.filter(pk__lt=obj.pk).order_by("-pk")
        if previousobj:
            return previousobj[0]
        else:
            return None

    def getPreviousObjEntries(self, queryset: list, obj: DjangoObjectType):
        """Finds all previous objects relative to obj of queryset"""
        previousqueryset = queryset.filter(pk__lt=obj.pk).order_by("-pk")
        return previousqueryset

    def getNextObjEntries(self, queryset: list, obj: DjangoObjectType):
        """Finds all next objects relative to obj of queryset"""
        nextqueryset = queryset.filter(pk__gt=obj.pk).order_by("pk")
        return nextqueryset

    def checkPreviousReactionProducts(self, reaction_id: int, smiles: str):
        """Checks if any previous reactions had the product matching the smiles"""
        reactionobj = self.getReaction(reaction_id=reaction_id)
        reactionqueryset = self.getReactionQuerySet(method_id=reactionobj.method_id.id)
        prevreactionqueryset = self.getPreviousObjEntries(
            queryset=reactionqueryset, obj=reactionobj
        )
        productmatches = []
        if prevreactionqueryset:
            for reactionobj in prevreactionqueryset:
                productobj = self.getProduct(reaction_id=reactionobj)
                if productobj.smiles == smiles:
                    productmatches.append(productobj)
            if productmatches:
                return True
            else:
                return False
        else:
            return False

    def checkNextReactionsAddActions(self, reactionobj: DjangoObjectType, smiles: str):
        """Checks if there are any reaction obj following the reaction.
        This could be the very next next or multiple steps later. Could possibly be used for
        several steps as well!
        If there is, gets AddAction matching smiles
        Args:
             reaction_id (int): Reaction id to use to search for next Reactions
             smiles (str): Product smiles of previous reaction to match with AddAction of
                           next Reaction
        """
        reactionqueryset = self.getReactionQuerySet(method_id=reactionobj.method_id.id)
        nextreactionqueryset = self.getNextObjEntries(
            queryset=reactionqueryset, obj=reactionobj
        )
        addactionsmatches = []
        for reactionobj in nextreactionqueryset:
            addactionmatch = self.getAddActions(reaction_id=reactionobj.id).filter(
                materialsmiles=smiles
            )
            if addactionmatch:
                addactionsmatches.append(addactionmatch[0])
        return addactionsmatches

    def getAnalyseActionQuerySet(self, reaction_ids: list):
        """Get Analyse Action queryset for list reaction ids object matching reaction_id and smiles
        Args:
            reaction_ids (list): List Reaction ids to search for
        Returns:
            analyseactionqueryset (Django_objs): AnalyseAction queryset
        """
        analyseactionqueryset = AnalyseAction.objects.filter(
            reaction_id__in=reaction_ids,
        )
        return analyseactionqueryset

    def getReaction(self, reaction_id: int):
        """Get reaction object
        Args:
            reaction_id (int): Reaction DB id to search for
        Returns:
            reactionobj (Django_obj): Reaction Django object
        """
        reactionobj = Reaction.objects.get(id=reaction_id)
        return reactionobj

    def getReactionQuerySet(self, method_id: int):
        """Get reaction queryset for method_id
        Args:
            method_id (int): Method DB id to search for
        Returns:
            reactionqueryset (Django_queryset): Reaction queryset related to method_id
        """
        reactionqueryset = Reaction.objects.filter(method_id=method_id)
        return reactionqueryset

    def getProductQuerySet(self, reaction_ids: list):
        """Get Product queryset from list of reaction_ids
        Args:
            reaction_ids (list): Reaction DB ids to search for
        Returns:
            productqueryset (Django_obj): Product queryset related to reaction_ids
        """
        productqueryset = Product.objects.filter(reaction_id__in=reaction_ids)
        return productqueryset

    def getProduct(self, reaction_id: int):
        """Get product object
        Args:
            reaction_id (int): Reaction DB id to search for
        Returns:
            productobj (Django_obj): Product Django object related to reaction_id
        """
        productobj = Product.objects.get(reaction_id=reaction_id)
        return productobj

    def getAddActions(self, reaction_id: int):
        """Get add actions queryset for reaction_id
        Args:
            reaction_id (int): Reaction DB id to search for
        Returns:
            addactionqueryset (Django_queryset): Reaction Django object
        """
        addactionqueryset = AddAction.objects.filter(reaction_id=reaction_id).order_by(
            "id"
        )
        return addactionqueryset

    def getStartingReactionPlateQuerySet(self):
        reactionplatequeryset = Plate.objects.filter(
            otsession_id=self.otsessionobj.id, name__contains="Reactionplate"
        ).order_by("id")
        return reactionplatequeryset

    def getRoundedAnalyseActionVolumes(self, analyseactionqueryset):
        roundedvolumes = [
            round(analyseactionobj.samplevolume + analyseactionobj.solventvolume)
            for analyseactionobj in analyseactionqueryset
        ]
        return roundedvolumes

    def getRoundedAddActionVolumes(self, addactionqueryset):
        roundedvolumes = [
            round(addactionobj.volume) for addactionobj in addactionqueryset
        ]
        return roundedvolumes

    def getTipRackType(self, roundedvolumes):
        tipsavailable = {
            300: "opentrons_96_tiprack_300ul",
            10: "opentrons_96_tiprack_20ul",
        }
        tipkey = min(
            tipsavailable,
            key=lambda x: self.getNumberTransfers(
                pipettevolume=x, roundedvolumes=roundedvolumes
            ),
        )
        tipracktype = tipsavailable[tipkey]
        return tipracktype

    def getNumberTips(self, queryset):
        numbertips = len(queryset)
        return numbertips

    def getNumberTransfers(self, pipettevolume, roundedvolumes):
        numbertransfers = sum(
            [
                round(volume / pipettevolume) if pipettevolume < volume else 1
                for volume in roundedvolumes
            ]
        )
        return numbertransfers

    def getPipetteType(self, roundedvolumes):
        pipettesavailable = {
            10: {
                "labware": "p10_single",
                "position": "right",
                "type": "single",
                "maxvolume": 10,
            },
            300: {
                "labware": "p300_single",
                "position": "right",
                "type": "single",
                "maxvolume": 300,
            },
        }
        pipettekey = min(
            pipettesavailable,
            key=lambda x: self.getNumberTransfers(
                pipettevolume=x, roundedvolumes=roundedvolumes
            ),
        )
        pipettetype = pipettesavailable[pipettekey]
        return pipettetype

    def getAnalyseActionsDataFrame(self):
        # Optimise -> https://stackoverflow.com/questions/11697887/converting-django-queryset-to-pandas-dataframe
        analyseactionsdf = pd.DataFrame(list(self.allanalyseactionqueryset.values()))
        analyseactionsdf = analyseactionsdf.rename(columns={"solventvolume": "volume"})
        return analyseactionsdf

    def getAddActionsDataFrame(self):
        # Optimise -> https://stackoverflow.com/questions/11697887/converting-django-queryset-to-pandas-dataframe
        addactionslistdf = []

        for addactionqueryset in self.alladdactionqueryset:
            addactionsdf = pd.DataFrame(list(addactionqueryset.values()))
            if not addactionsdf.empty:
                addactionslistdf.append(addactionsdf)
        addactionsdf = pd.concat(addactionslistdf)

        addactionsdf["uniquesolution"] = addactionsdf.apply(
            lambda row: self.combinestrings(row), axis=1
        )

        return addactionsdf

    def getMaxWellVolume(self, plateobj):
        maxwellvolume = plateobj.maxwellvolume
        return maxwellvolume

    def getDeadVolume(self, maxwellvolume):
        deadvolume = maxwellvolume * 0.05
        return deadvolume

    def getCloneWells(self, plateobj):
        clonewellqueryset = Well.objects.filter(plate_id=plateobj.id)
        return clonewellqueryset

    def getUniqueAnalyseMethods(self):
        methods = sorted(
            set(
                [
                    analyseactionobj.method
                    for analyseactionobj in self.allanalyseactionqueryset
                ]
            )
        )
        return methods

    def getUniqueTemperatures(self):
        temperatures = sorted(
            set([reactionobj.temperature for reactionobj in self.reactiongroupqueryset])
        )
        return temperatures

    # def getDistinctFieldValues(self, djangomodel: Django_model, field: str):
    #     """Returns distinct values of a a field for a Django model
    #     """
    #     distinctvaluesqueryset = djangomodel.__class__.objects.order_by().values_list(field, flat=True).distinct()
    #     return distinctvaluesqueryset

    # def getGroupedObjects(self, queryset: list, fieldvalues: list, field: str):
    #     """Group objects by field value of a input queryset and fieldvalues
    #     """
    #     # djangomodel = queryset.model
    #     # distinctvaluesqueryset = self.getDistinctFieldValues(djangomodel=djangomodel, field=field)
    #     groupedobjs = []

    #     for fieldvalue in fieldvalues:
    #         group = [
    #             obj
    #             for obj in queryset
    #             if getattr(obj, field) == fieldvalue
    #         ]
    #         groupedobjs.append(group)

    #     return groupedobjs

    def getGroupedAnalyseMethods(self):
        """Group Analyse methods"""
        methods = self.getUniqueAnalyseMethods()
        groupedanalysemethodobjs = []

        for method in methods:
            analysemethodgroup = [
                analyseactionobj
                for analyseactionobj in self.allanalyseactionqueryset
                if analyseactionobj.method == method
            ]
            groupedanalysemethodobjs.append(analysemethodgroup)

        return groupedanalysemethodobjs

    def getGroupedTemperatureReactions(self):
        """Group reactions done at the same temperature"""
        temperatures = self.getUniqueTemperatures()
        groupedtemperaturereactionobjs = []

        for temperature in temperatures:
            temperaturereactiongroup = [
                reactionobj
                for reactionobj in self.reactiongroupqueryset
                if reactionobj.temperature == temperature
            ]
            groupedtemperaturereactionobjs.append(temperaturereactiongroup)

        return groupedtemperaturereactionobjs

    def getMedianValue(self, values):
        medianvalue = median(values)
        return medianvalue

    def getMaxValue(self, values):
        maxvalue = max(values)
        return maxvalue

    def getSumValue(self, values):
        sumvalue = sum(values)
        return sumvalue

    def getNumberObjs(self, queryset: list):
        numberobjs = len(queryset)
        return numberobjs

    def getReactionLabwarePlateType(self, grouptemperaturereactionobjs):
        numberreactions = self.getNumberObjs(queryset=grouptemperaturereactionobjs)
        reactionvolumes = []
        for reactionobj in grouptemperaturereactionobjs:
            addactionqueryset = self.getAddActions(reaction_id=reactionobj.id)
            roundedvolumes = self.getRoundedAddActionVolumes(
                addactionqueryset=addactionqueryset
            )
            sumvolume = self.getSumValue(values=roundedvolumes)
            reactionvolumes.append(sumvolume)
        maxvolume = self.getMaxValue(values=reactionvolumes)
        medianvolume = self.getMedianValue(values=reactionvolumes)
        headspacevolume = maxvolume + (maxvolume * 0.2)
        reactiontemperature = grouptemperaturereactionobjs[0].temperature

        if reactiontemperature > 25:
            labwareplatetypes = [
                labware_plate
                for labware_plate in labware_plates
                if labware_plates[labware_plate]["reflux"]
                and labware_plates[labware_plate]["volume_well"] > headspacevolume
                and labware_plates[labware_plate]["no_wells"] >= numberreactions
            ]
        else:
            labwareplatetypes = [
                labware_plate
                for labware_plate in labware_plates
                if not labware_plates[labware_plate]["reflux"]
                and labware_plates[labware_plate]["volume_well"] > headspacevolume
                and labware_plates[labware_plate]["no_wells"] >= numberreactions
            ]

        if len(labwareplatetypes) > 1:
            volumewells = [
                labware_plates[labware_plate]["volume_well"]
                for labware_plate in labwareplatetypes
            ]
            indexclosestvalue = min(
                range(len(volumewells)), key=lambda x: abs(x - medianvolume)
            )
            labwareplatetype = labwareplatetypes[indexclosestvalue]
        else:
            labwareplatetype = labwareplatetypes[0]

        return labwareplatetype

    def getNumberVials(self, maxvolumevial, volumematerial):
        if maxvolumevial > volumematerial:
            novialsneeded = 1
        else:
            volumestoadd = []
            deadvolume = self.getDeadVolume(maxwellvolume=maxvolumevial)
            novialsneededratio = volumematerial / (maxvolumevial - deadvolume)
            frac, whole = math.modf(novialsneededratio)
            volumestoadd = [maxvolumevial for maxvolumevial in range(int(whole))]
            volumestoadd.append(frac * maxvolumevial + deadvolume)
            novialsneeded = sum(volumestoadd)
        return novialsneeded

    def getAnalysePlateType(self, wellsneeded: int, method: str):
        """Finds best labware plate for analysis"""
        analyselc = method.lower()

        labwareplatetypes = [
            labware_plate
            for labware_plate in labware_plates
            if labware_plates[labware_plate]["{}plate".format(analyselc)]
        ]

        wellcomparedict = {}

        for labwareplate in labwareplatetypes:
            maxvolumevial = labware_plates[labwareplate]["volume_well"]
            noplatevials = labware_plates[labwareplate]["no_wells"]
            platesneeded = int(math.ceil(wellsneeded / noplatevials))
            wellcomparedict[maxvolumevial] = platesneeded
            minimumnovialsvolume = min(wellcomparedict, key=wellcomparedict.get)

        labwareplatetype = [
            labware_plate
            for labware_plate in labwareplatetypes
            if labware_plates[labware_plate]["volume_well"] == minimumnovialsvolume
        ][0]

        return labwareplatetype

    def getStarterPlateType(self, startingmaterialsdf: DataFrame):
        labwareplatetypes = [
            labware_plate
            for labware_plate in labware_plates
            if labware_plates[labware_plate]["starting_plate"]
        ]

        vialcomparedict = {}

        for labwareplate in labwareplatetypes:
            maxvolumevial = labware_plates[labwareplate]["volume_well"]
            noplatevials = labware_plates[labwareplate]["no_wells"]

            vialsneeded = startingmaterialsdf.apply(
                lambda row: self.getNumberVials(
                    maxvolumevial=maxvolumevial, volumematerial=row["volume"]
                ),
                axis=1,
            )
            totalvialsneeded = sum(vialsneeded)
            platesneeded = int(math.ceil(totalvialsneeded / noplatevials))

            vialcomparedict[maxvolumevial] = platesneeded

        minimumnovialsvolume = min(vialcomparedict, key=vialcomparedict.get)

        labwareplatetype = [
            labware_plate
            for labware_plate in labwareplatetypes
            if labware_plates[labware_plate]["volume_well"] == minimumnovialsvolume
        ][0]

        return labwareplatetype

    def getAddActionsMaterialDataFrame(self, productexists: bool):
        """Aggegates all materials and sums up volume requires using solvent type and
        concentration
        Args:
             productexists (bool): Set to true to get materials that are products
                                   of the previous reaction step
        """
        materialsdf = self.addactionsdf.groupby(["uniquesolution"]).agg(
            {
                "reaction_id_id": "first",
                "smiles": "first",
                "volume": "sum",
                "solvent": "first",
                "concentration": "first",
            }
        )
        materialsdf["productexists"] = materialsdf.apply(
            lambda row: self.checkPreviousReactionProducts(
                reaction_id=row["reaction_id_id"], smiles=row["smiles"]
            ),
            axis=1,
        )

        if productexists:
            materialsdf = materialsdf[materialsdf["productexists"]]

        if not productexists:
            materialsdf = materialsdf[~materialsdf["productexists"]]

        materialsdf = materialsdf.sort_values(["solvent", "volume"], ascending=False)

        return materialsdf

    def createOTSessionModel(self):
        otsessionobj = OTSession()
        otsessionobj.otbatchprotocol_id = self.otbatchprotocolobj
        otsessionobj.reactionstep = self.reactionstep
        otsessionobj.sessiontype = self.sessiontype
        otsessionobj.save()
        return otsessionobj

    def createDeckModel(self):
        deckobj = Deck()
        deckobj.otsession_id = self.otsessionobj
        deckobj.numberslots = 11
        deckobj.save()
        self.deckobj = deckobj
        return deckobj

    def createPipetteModel(self):
        pipetteobj = Pipette()
        pipetteobj.otsession_id = self.otsessionobj
        pipetteobj.position = self.pipettetype["position"]
        pipetteobj.maxvolume = self.pipettetype["maxvolume"]
        pipetteobj.type = self.pipettetype["type"]
        pipetteobj.name = "{}_{}".format(
            self.pipettetype["position"], self.pipettetype["labware"]
        )
        pipetteobj.labware = self.pipettetype["labware"]
        pipetteobj.save()

    def createTiprackModel(self, name: str):
        """Creates TipRack model"""
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            index = indexslot
            tiprackobj = TipRack()
            tiprackobj.otsession_id = self.otsessionobj
            tiprackobj.deck_id = self.deckobj
            tiprackobj.name = "{}_{}".format(name, indexslot)
            tiprackobj.index = index
            tiprackobj.labware = name
            tiprackobj.save()
        else:
            print("No more deck slots available")

    def createPlateModel(self, platetype: str, platename: str, labwaretype: str):
        """Creates Plate model if Deck index is available"""
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            plateindex = indexslot
            maxwellvolume = labware_plates[labwaretype]["volume_well"]
            numberwells = labware_plates[labwaretype]["no_wells"]
            plateobj = Plate()
            plateobj.otsession_id = self.otsessionobj
            plateobj.deck_id = self.deckobj
            plateobj.labware = labwaretype
            plateobj.index = plateindex
            plateobj.name = "Reaction_step_{}_{}_index_{}".format(
                self.reactionstep, platename, indexslot
            )
            plateobj.type = platetype
            plateobj.maxwellvolume = maxwellvolume
            plateobj.numberwells = numberwells
            plateobj.save()
            return plateobj
        else:
            print("No more deck slots available")

    def createWellModel(
        self,
        plateobj,
        welltype,
        wellindex,
        volume=None,
        reactionobj=None,
        smiles=None,
        concentration=None,
        solvent=None,
        reactantfornextstep=False,
    ):
        wellobj = Well()
        wellobj.otsession_id = self.otsessionobj
        wellobj.plate_id = plateobj
        if reactionobj:
            wellobj.reaction_id = reactionobj
            wellobj.method_id = reactionobj.method_id
        wellobj.type = welltype
        wellobj.index = wellindex
        wellobj.volume = volume
        wellobj.smiles = smiles
        wellobj.concentration = concentration
        wellobj.solvent = solvent
        wellobj.reactantfornextstep = reactantfornextstep
        wellobj.save()
        return wellobj

    def calcMass(self, row):
        mols = row["concentration"] * row["amount-ul"] * 1e-6
        smiles = row["SMILES"]
        mw = Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles))
        massmg = mols * mw * 1e3
        return round(massmg, 2)

    def createCompoundOrderModel(self, orderdf):
        compoundorderobj = CompoundOrder()
        compoundorderobj.otsession_id = self.otsessionobj
        csvdata = orderdf.to_csv(encoding="utf-8", index=False)
        ordercsv = default_storage.save(
            "compoundorders/"
            + "{}-session-starterplate-for-batch-{}-reactionstep-{}-sessionid-{}".format(
                self.sessiontype,
                self.batchobj.batchtag,
                self.reactionstep,
                str(self.otsessionobj.id),
            )
            + ".csv",
            ContentFile(csvdata),
        )
        compoundorderobj.ordercsv = ordercsv
        compoundorderobj.save()

    def createSolventPrepModel(self, solventdf):
        solventprepobj = SolventPrep()
        solventprepobj.otsession_id = self.otsessionobj
        csvdata = solventdf.to_csv(encoding="utf-8", index=False)
        ordercsv = default_storage.save(
            "solventprep/"
            + "{}-session-solventplate-for-batch-{}-reactionstep-{}-sessionid-{}".format(
                self.sessiontype,
                self.batchobj.batchtag,
                self.reactionstep,
                str(self.otsessionobj.id),
            )
            + ".csv",
            ContentFile(csvdata),
        )
        solventprepobj.solventprepcsv = ordercsv
        solventprepobj.save()

    def createTipRacks(self, tipracktype: str):
        numberacks = int(-(-self.numbertips // 96))
        for rack in range(numberacks):
            self.createTiprackModel(name=tipracktype)

    def checkDeckSlotAvailable(self):
        testslotavailable = self.deckobj.indexslotavailable
        if testslotavailable <= self.deckobj.numberslots:
            self.deckobj.indexslotavailable = testslotavailable + 1
            self.deckobj.save()
            return testslotavailable
        else:
            self.deckobj.slotavailable = False
            self.deckobj.save()
            return False

    def checkPlateWellsAvailable(self, plateobj):
        wellavailable = plateobj.indexswellavailable + 1
        numberwells = plateobj.numberwells
        if wellavailable <= numberwells:
            plateobj.indexswellavailable = wellavailable
            plateobj.save()
            return plateobj.indexswellavailable
        else:
            plateobj.wellavailable = False
            plateobj.save()
            return False

    def createAnalysePlate(
        self, method: str, platetype: str, analyseactionqueryset: list
    ):
        plateobj = self.createPlateModel(
            platetype="analyse",
            platename="{}_analyse_plate".format(method),
            labwaretype=platetype,
        )

        for analyseactionobj in analyseactionqueryset:
            reactionobj = analyseactionobj.reaction_id
            productobj = self.getProduct(reaction_id=reactionobj.id)
            volume = analyseactionobj.samplevolume + analyseactionobj.solventvolume
            indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
            if not indexwellavailable:
                plateobj = self.createPlateModel(
                    platename="{}_analyse_plate".format(method),
                    labwaretype=platetype,
                )

                indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)

            self.createWellModel(
                plateobj=plateobj,
                reactionobj=reactionobj,
                welltype="analyse",
                wellindex=indexwellavailable - 1,
                volume=volume,
                smiles=productobj.smiles,
            )

    def createReactionStartingPlate(self):
        startingmaterialsdf = self.getAddActionsMaterialDataFrame(productexists=False)

        startinglabwareplatetype = self.getStarterPlateType(
            startingmaterialsdf=startingmaterialsdf
        )

        plateobj = self.createPlateModel(
            platetype="startingmaterial",
            platename="Startingplate",
            labwaretype=startinglabwareplatetype,
        )

        maxwellvolume = self.getMaxWellVolume(plateobj=plateobj)
        deadvolume = self.getDeadVolume(maxwellvolume=maxwellvolume)

        orderdictslist = []

        for i in startingmaterialsdf.index.values:
            extraerrorvolume = startingmaterialsdf.at[i, "volume"] * 0.1
            totalvolume = startingmaterialsdf.at[i, "volume"] + extraerrorvolume
            if totalvolume > maxwellvolume:
                nowellsneededratio = totalvolume / (maxwellvolume - deadvolume)

                frac, whole = math.modf(nowellsneededratio)
                volumestoadd = [maxwellvolume for i in range(int(whole))]
                volumestoadd.append(frac * maxwellvolume + deadvolume)

                for volumetoadd in volumestoadd:
                    indexwellavailable = self.checkPlateWellsAvailable(
                        plateobj=plateobj
                    )
                    if not indexwellavailable:

                        plateobj = self.createPlateModel(
                            platetype="startingmaterial",
                            platename="Startingplate",
                            labwaretype=startinglabwareplatetype,
                        )

                        indexwellavailable = self.checkPlateWellsAvailable(
                            plateobj=plateobj
                        )

                    wellobj = self.createWellModel(
                        plateobj=plateobj,
                        reactionobj=self.getReaction(
                            reactionid=startingmaterialsdf.at[i, "reaction_id_id"]
                        ),
                        welltype="startingmaterial",
                        wellindex=indexwellavailable - 1,
                        volume=volumetoadd,
                        smiles=startingmaterialsdf.at[i, "smiles"],
                        concentration=startingmaterialsdf.at[i, "concentration"],
                        solvent=startingmaterialsdf.at[i, "solvent"],
                    )

                    orderdictslist.append(
                        {
                            "SMILES": startingmaterialsdf.at[i, "smiles"],
                            "platename": plateobj.name,
                            "well": wellobj.index,
                            "concentration": startingmaterialsdf.at[i, "concentration"],
                            "solvent": startingmaterialsdf.at[i, "solvent"],
                            "amount-ul": round(volumetoadd, 2),
                        }
                    )

            else:
                indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                volumetoadd = totalvolume + deadvolume

                if not indexwellavailable:
                    plateobj = self.createPlateModel(
                        platetype="startingmaterial",
                        platename="Startingplate",
                        labwaretype=startinglabwareplatetype,
                    )
                    indexwellavailable = self.checkPlateWellsAvailable(
                        plateobj=plateobj
                    )

                wellobj = self.createWellModel(
                    plateobj=plateobj,
                    reactionobj=self.getReaction(
                        reaction_id=startingmaterialsdf.at[i, "reaction_id_id"]
                    ),
                    welltype="startingmaterial",
                    wellindex=indexwellavailable - 1,
                    volume=volumetoadd,
                    smiles=startingmaterialsdf.at[i, "smiles"],
                    concentration=startingmaterialsdf.at[i, "concentration"],
                    solvent=startingmaterialsdf.at[i, "solvent"],
                )

                orderdictslist.append(
                    {
                        "SMILES": startingmaterialsdf.at[i, "smiles"],
                        "platename": plateobj.name,
                        "well": wellobj.index,
                        "concentration": startingmaterialsdf.at[i, "concentration"],
                        "solvent": startingmaterialsdf.at[i, "solvent"],
                        "amount-ul": round(volumetoadd, 2),
                    }
                )

        orderdf = pd.DataFrame(orderdictslist)
        orderdf["mass-mg"] = orderdf.apply(lambda row: self.calcMass(row), axis=1)

        self.createCompoundOrderModel(orderdf=orderdf)

    def createReactionPlate(self):
        for grouptemperaturereactionobjs in self.groupedtemperaturereactionobjs:
            labwareplatetype = self.getReactionLabwarePlateType(
                grouptemperaturereactionobjs=grouptemperaturereactionobjs
            )

            plateobj = self.createPlateModel(
                platetype="reaction",
                platename="Reactionplate",
                labwaretype=labwareplatetype,
            )

            for reactionobj in grouptemperaturereactionobjs:
                productobj = self.getProduct(reaction_id=reactionobj.id)
                indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                if not indexwellavailable:
                    plateobj = self.createPlateModel(
                        platetype="reaction",
                        platename="Reactionplate",
                        labwaretype=labwareplatetype,
                    )

                    indexwellavailable = self.checkPlateWellsAvailable(
                        plateobj=plateobj
                    )

                nextaddactionobjs = self.checkNextReactionsAddActions(
                    reactionobj=reactionobj, smiles=productobj.smiles
                )
                if nextaddactionobjs:
                    nextaddactionobj = nextaddactionobjs[0]
                    self.createWellModel(
                        plateobj=plateobj,
                        reactionobj=reactionobj,
                        welltype="reaction",
                        wellindex=indexwellavailable - 1,
                        volume=nextaddactionobj.volume,
                        smiles=productobj.smiles,
                        concentration=nextaddactionobj.concentration,
                        solvent=nextaddactionobj.solvent,
                        reactantfornextstep=True,
                    )

                else:
                    self.createWellModel(
                        plateobj=plateobj,
                        reactionobj=reactionobj,
                        welltype="reaction",
                        wellindex=indexwellavailable - 1,
                        smiles=productobj.smiles,
                    )

    def combinestrings(self, row):
        return (
            str(row["smiles"])
            + "-"
            + str(row["solvent"])
            + "-"
            + str(row["concentration"])
        )

    def cloneInputPlate(self):
        for plateobj in self.inputplatequeryset:
            indexslot = self.checkDeckSlotAvailable()
            if indexslot:
                clonewellqueryset = self.getCloneWells(plateobj=plateobj)
                plateindex = indexslot
                previousname = plateobj.name
                platename = "Startingplate"
                plateobj.pk = None
                plateobj.deck_id = self.deckobj
                plateobj.otsession_id = self.otsessionobj
                plateobj.index = plateindex
                plateobj.name = "{}_{}_from_{}".format(
                    platename, indexslot, previousname
                )
                plateobj.save()
                self.cloneInputWells(clonewellqueryset, plateobj)
            else:
                print("No more deck slots available")

    def cloneInputWells(self, clonewellqueryset, plateobj):
        for clonewellobj in clonewellqueryset:
            clonewellobj.pk = None
            clonewellobj.plate_id = plateobj
            clonewellobj.otsession_id = self.otsessionobj
            clonewellobj.save()

    def createSolventPlate(self, materialsdf: DataFrame):
        """Creates solvent plate/s for diluting reactants for reactions or analysis."""
        if not materialsdf.empty:
            solventdictslist = []
            materialsdf = materialsdf.groupby(["solvent"])["volume"].sum().to_frame()
            startinglabwareplatetype = self.getStarterPlateType(
                startingmaterialsdf=materialsdf
            )

            plateobj = self.createPlateModel(
                platetype="dilution",
                platename="Solventplate",
                labwaretype=startinglabwareplatetype,
            )

            maxwellvolume = self.getMaxWellVolume(plateobj=plateobj)
            deadvolume = self.getDeadVolume(maxwellvolume=maxwellvolume)

            for solventgroup in materialsdf.index.values:
                totalvolume = materialsdf.at[solventgroup, "volume"]
                if totalvolume > maxwellvolume:
                    nowellsneededratio = totalvolume / (maxwellvolume - deadvolume)
                    frac, whole = math.modf(nowellsneededratio)
                    volumestoadd = [maxwellvolume for i in range(int(whole))]
                    volumestoadd.append(frac * maxwellvolume + deadvolume)

                    for volumetoadd in volumestoadd:
                        indexwellavailable = self.checkPlateWellsAvailable(
                            plateobj=plateobj
                        )
                        if not indexwellavailable:

                            plateobj = self.createPlateModel(
                                platetype="dilution",
                                platename="Solventplate",
                                labwaretype=startinglabwareplatetype,
                            )

                            indexwellavailable = self.checkPlateWellsAvailable(
                                plateobj=plateobj
                            )

                        wellobj = self.createWellModel(
                            plateobj=plateobj,
                            welltype="dilution",
                            wellindex=indexwellavailable - 1,
                            volume=volumetoadd,
                            solvent=solventgroup,
                        )

                        solventdictslist.append(
                            {
                                "platename": plateobj.name,
                                "well": wellobj.index,
                                "solvent": solventgroup,
                                "amount-ul": volumetoadd,
                            }
                        )

                else:
                    indexwellavailable = self.checkPlateWellsAvailable(
                        plateobj=plateobj
                    )
                    volumetoadd = totalvolume + deadvolume

                    if not indexwellavailable:
                        plateobj = self.createPlateModel(
                            platetype="dilution",
                            platename="Solventplate",
                            labwaretype=startinglabwareplatetype,
                        )
                        indexwellavailable = self.checkPlateWellsAvailable(
                            plateobj=plateobj
                        )

                    wellobj = self.createWellModel(
                        plateobj=plateobj,
                        welltype="dilution",
                        wellindex=indexwellavailable - 1,
                        volume=volumetoadd,
                        solvent=solventgroup,
                    )

                    solventdictslist.append(
                        {
                            "platename": plateobj.name,
                            "well": wellobj.index,
                            "solvent": solventgroup,
                            "amount-ul": volumetoadd,
                        }
                    )

            solventdf = pd.DataFrame(solventdictslist)

            self.createSolventPrepModel(solventdf=solventdf)
