"""Create LCMS starter csv for a batch"""
import math
from graphene_django import DjangoObjectType
from car.models import Deck, LCMSSession, OTBatchProtocol, Plate, Well
import pandas as pd


class CreateLCMSSession(object):
    """ "
    Creates an input csv for LCMS analysis for OT session
    """

    def __init__(
        self,
        protocolname: str,
        otsessionobj: DjangoObjectType,
        projectobj: DjangoObjectType,
        maximumlcmswells: int,
        maximumlcmsplates: int,
    ):
        self.protocolname = protocolname
        self.otsessionobj = otsessionobj
        self.projectobj = projectobj
        self.maximumlcmswells = maximumlcmswells
        self.maximwellsperplate = maximumlcmswells / maximumlcmsplates
        self.reactionplatequeryset = self.getReactionPlateQuerySet()
        self.reactionplateids = [
            reactionplate.id for reactionplate in self.reactionplatequeryset
        ]
        self.wellqueryset = self.getWellQuerySet()
        if self.checkMaxNoWells():
            self.wellstobeanalysed = self.wellqueryset[0 : self.maximumlcmswells]
            self.wellsfornextlcmssession = self.wellqueryset[self.maximumlcmswells : :]
            self.lcmssessionobj = self.createLCMSSession()
            self.createDeckModel()

        if not self.checkMaxNoWells():
            self.wellstobeanalysed = self.wellqueryset
            self.lcmssessionobj = self.createLCMSSession()
            self.createDeckModel()

    def getNoPlates(self):
        totalwellneeded = len(self.wellqueryset)
        platesneeded = int(math.ceil(totalwellneeded / "XX"))

    def createLCMSSession(self):
        lcmssession = LCMSSession()
        lcmssession.project_id = self.projectobj
        lcmssession.otsession_id = self.otsessionobj
        lcmssession.save()
        return lcmssession

    def createDeckModel(self):
        deckobj = Deck()
        deckobj.otsession_id = self.otsessionobj
        deckobj.lcmssession_id = self.lcmssessionobj
        deckobj.numberslots = 2
        deckobj.save()
        self.deckobj = deckobj
        return deckobj

    def createPlateModel(self, platename, labwaretype):
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            plateindex = indexslot
            maxwellvolume = labware_plates[labwaretype]["volume_well"]
            numberwells = labware_plates[labwaretype]["no_wells"]
            plateobj = Plate()
            plateobj.otsession_id = self.otsessionobj
            plateobj.deck_id = self.deckobj
            plateobj.platename = "Reaction_step_{}_{}_index_{}".format(
                self.reactionstep, platename, indexslot
            )
            plateobj.plateindex = plateindex
            plateobj.labware = labwaretype
            plateobj.maxwellvolume = maxwellvolume
            plateobj.numberwells = numberwells
            plateobj.save()
            return plateobj
        else:
            print("No more deck slots available")

    def createWellModel(
        self,
        plateobj,
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
        wellobj.wellindex = wellindex
        wellobj.volume = volume
        wellobj.smiles = smiles
        wellobj.concentration = concentration
        wellobj.solvent = solvent
        wellobj.reactantfornextstep = reactantfornextstep
        wellobj.save()
        return wellobj

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

    def getReactionPlateQuerySet(self):
        reactionplatequeryset = Plate.objects.filter(
            otsession_id=self.otsessionobj.id, platename="Reactionplate"
        )
        return reactionplatequeryset

    def getWellQuerySet(self, plate_ids: int):
        """Get Well queryset using list of reaction plate_ids"""
        wellqueryset = Well.objects.filter(plate_id__in=self.reactionplateids)
        return wellqueryset

    def checkMaxNoWells(self):
        """Checks if the number of wells needed for analysis is more than
        2 x 384 plates. The max number of plates allowed on Agilent
        LCMS is 2.
        """
        if len(self.wellqueryset) > self.maximumlcmswells:
            return True
        else:
            return False

    def createLCMSdf(self):
        lcmsdf = pd.DataFrame(
            columns=[
                "Line",
                "Location",
                "Sample Name",
                "Method Name",
                "Inj/Location",
                "Sample Type",
                "inj Volume",
                "Datafile",
            ]
        )
        return lcmsdf

    def createLCMSPlate(self):
        if len(self.wellstobeanalysed) < self.maximwellsperplate:
            lcmsdf = self.createLCMSdf()
            lcmsdf["Line"] = [
                index + 1 for index, value in enumerate(self.wellstobeanalysed)
            ]
            lcmsdf["Sample Name"] = [well.smiles for well in self.wellstobeanalysed]
            lcmsdf["Method Name"] = self.lcmsmethodname

    def createFilePath(self):
        filename = "LCMS-inputcsv-batch-{}-reactionstep{}-sessionid-{}.txt".format(
            self.protocolname, self.reactionstep, self.otsessionid
        )
        path = "tmp/" + filename
        filepath = str(os.path.join(settings.MEDIA_ROOT, path))
        return filepath, filename
