from pyexpat import model
from django.db import models


class Project(models.Model):
    """Django model to define a Project - a compound synthesis project.

    Parameters
    ----------
    init_date: DateTimeField
        The date the target was initiated (autofield)
    name: SlugFieldField
        The name of the project created as combination of the first three letters of the
        submitter's name and submitter's organisation
    submittername: Charfield
        The name of the person creating the project
    submitterorganisation: Charfield
        The name of the person's organisation creating the project
    proteintarget: Charfield
        The target protein for the submitted compounds
    quotecost: FloatField
        Optional field for the total cost of the project using MCule as a supplier
    quoteurl: CharField
        Optional field for the MCule quote url
    """

    init_date = models.DateTimeField(auto_now_add=True)
    name = models.SlugField(max_length=100, db_index=True)
    submitterorganisation = models.CharField(max_length=100)
    submittername = models.CharField(max_length=255)
    proteintarget = models.CharField(max_length=100)
    quotedcost = models.FloatField(null=True)
    quoteurl = models.CharField(max_length=255, null=True)


class Batch(models.Model):
    """Django model to define a Batch - a batch of compounds.

    Parameters
    ----------
    project_id: ForeignKey
        Foreign key linking a batch to it's project
    batch_id: ForeignKey
        Optional foreign key linking a sub-batch to it's parent-batch
    batchtag: Charfield
        The name of the batch
    """

    project_id = models.ForeignKey(
        Project, related_name="batches", on_delete=models.CASCADE
    )
    batch_id = models.ForeignKey("Batch", on_delete=models.CASCADE, null=True)
    batchtag = models.CharField(max_length=50)


class Target(models.Model):
    """Django model to define a Target - a target compound for synthesis.

    Parameters
    ----------
    batch_id: ForeignKey
        Foreign key linking a target to it's batch
    smiles: Charfield
        The SMILES of the target compound
    image: FileField
        File link to a stored version of the image file of the target compound
    name: Charfield
        The name of the compound (Must complete how this is made)
    mass: FloatField
        The target synthesis mass (mg) for a target compound
    mols: FloatField
        The target mols of the compound to be synthesised
    """

    batch_id = models.ForeignKey(
        Batch, related_name="targets", on_delete=models.CASCADE
    )
    smiles = models.CharField(max_length=255, db_index=True)
    image = models.FileField(upload_to="targetimages/", max_length=255)
    name = models.CharField(max_length=255, db_index=True)
    mass = models.FloatField()
    mols = models.FloatField()


class Method(models.Model):
    """Django model to define a Method - a retrosynthetic pathway of reactions
    for a target compound.

    Parameters
    ----------
    target_id: ForeignKey
        Foreign key linking a method to it's target compound
    nosteps: IntegerField
        The number of reaction steps for a method
    otchem: BooleanField
        Set to True if all the reactions in the method can be executed on the OpenTrons
    """

    target_id = models.ForeignKey(
        Target, related_name="methods", on_delete=models.CASCADE
    )
    nosteps = models.IntegerField()
    otchem = models.BooleanField(default=False)


class Reaction(models.Model):
    """Django model to define a Reaction - the reaction to make a compound.

    Parameters
    ----------
    method_id: ForeignKey
        Foreign key linking a reaction to it's method
    reactionclass: Charfield
        The name of the reaction
    intramolecular: BooleanField
        Set to True if the reaction is intermolecular
    recipetype: CharField
        The encoded recipe type used to execute the reaction
    image: FileField
        File link to a stored version of the image file of the reaction
    success: BooleanField
        The success of the reaction with default success set to True
    """

    method_id = models.ForeignKey(
        Method, related_name="reactions", on_delete=models.CASCADE
    )
    reactionclass = models.CharField(max_length=255)
    intramolecular = models.BooleanField(default=False)
    recipetype = models.CharField(max_length=50, default="standard", null=True)
    temperature = models.IntegerField(default=25)
    image = models.FileField(
        upload_to="reactionimages/",
        max_length=255,
        null=True,
    )
    success = models.BooleanField(default=True)


class PubChemInfo(models.Model):
    """Django model to define PubChemInfo - the PubChem info for a compound.

    Parameters
    ----------
    compoundid: IntegerField
        The PubChem compound id of a compound
    summaryurl: Charfield
        The PubChem compound url
    lcssurl: Charfield
        The PubChem Laboratory Chemical Safety Summary url
    smiles: Charfield
        The SMILES of the compound
    cas: Charfield
        Optional CAS number (if found) for a compound
    """

    compoundid = models.IntegerField()
    summaryurl = models.CharField(max_length=255)
    lcssurl = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    cas = models.CharField(max_length=50, null=True)


class Reactant(models.Model):
    """Django model to define a Reactant - the reactant compound of a reaction.

    Parameters
    ----------
    reaction_id: ForeignKey
        Optional foreign key linking a reactant to it's reaction
    pubcheminfo_id: ForeignKey
        Optional foreign key linking a reactant to it's pubchem info (if found)
    smiles: Charfield
        The SMILES of the reactant compound
    """

    reaction_id = models.ForeignKey(
        Reaction, related_name="reactants", on_delete=models.CASCADE, null=True
    )
    pubcheminfo_id = models.ForeignKey(
        PubChemInfo,
        related_name="reactantpubcheminfo",
        on_delete=models.PROTECT,
        null=True,
    )
    smiles = models.CharField(max_length=255)


class Product(models.Model):
    """Django model to define a Product- the product compound of a reaction.

    Parameters
    ----------
    reaction_id: ForeignKey
        Foreign key linking a product to it's reaction
    pubcheminfo_id: ForeignKey
        Optional foreign key linking a product to it's pubchem info (if found)
    smiles: Charfield
        The SMILES of the product compound
    image: FileField
        File link to a stored version of the image file of the product compound
    """

    reaction_id = models.ForeignKey(
        Reaction, related_name="products", on_delete=models.CASCADE
    )
    pubcheminfo_id = models.ForeignKey(
        PubChemInfo,
        related_name="productpubcheminfo",
        on_delete=models.PROTECT,
        null=True,
    )
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to="productimages/", max_length=255)


class CatalogEntry(models.Model):
    """Django model to define a CatalogEntry- the catalog information for a compound.

    Parameters
    ----------
    reactant_id: ForeignKey
        Optional Foreign key linking a catalog entry to a reactant
    target_id: ForeignKey
        Optional Foreign key linking a catalog entry to a product
    vendor: Charfield
        The vendor/supplier of the compound
    catalogid: Charfield
        The catalog/vendor id of the compound
    priceinfo: Charfield
        The catalog price info ($) for a compound. This can be a range -> "$100 < 1k/g"
    upperprice: IntegerField
        Optional upper price ($) of a compound. Highest price if range given.
    leadtime: IntegerField
        Optional lead time (weeks) for a compound to be delivered
    """

    reactant_id = models.ForeignKey(
        Reactant, related_name="catalogentries", on_delete=models.CASCADE, null=True
    )
    target_id = models.ForeignKey(
        Target, related_name="catalogentries", on_delete=models.CASCADE, null=True
    )

    vendor = models.CharField(max_length=100)
    catalogid = models.CharField(max_length=50)
    priceinfo = models.CharField(max_length=50)
    upperprice = models.IntegerField(null=True)
    leadtime = models.IntegerField(null=True)


class ActionSession(models.Model):
    class Type(models.TextChoices):
        reaction = "reaction"
        stir = "stir"
        workup = "workup"
        analyse = "analyse"

    class Driver(models.TextChoices):
        human = "human"
        robot = "robot"

    reaction_id = models.ForeignKey(
        Reaction, related_name="actionsessions", on_delete=models.CASCADE
    )
    sessionnumber = models.IntegerField()
    type = models.CharField(choices=Type.choices, max_length=10)
    driver = models.CharField(
        choices=Driver.choices, default=Driver.robot, max_length=10
    )


class AddAction(models.Model):
    """Django model to define a AddAction - the add action details

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    fromplatetype: CharField
        The plate the add action is moving material from
    toplatetype: CharField
        The plate the add action is moviong material to
    smiles: CharField
        Optional SMILES of the material being added
    calcunit: CharField
        The unit used for calculating the amount to add eg. molar eq. (moleq)
        and mass eq. (masseq)
    volume: FloatField
        The volume being added
    volumeunit: CharField
        The unit of the volume being added (default=ul)
    molecularweight: FloatField
        The molecular weight of the compound being added
    solvent: CharField
        Optional solvent used to dilute the material being added
    concentration: FloatField
        Optional concentration of the material solution prepared
    """

    class CalcUnit(models.TextChoices):
        moleq = "moleq"
        masseq = "masseq"
        ul = "ul"

    class VolumeUnit(models.TextChoices):
        ul = "ul"
        ml = "ml"

    class PlateType(models.TextChoices):
        reaction = "reaction"
        workup1 = "workup1"
        workup2 = "workup2"
        workup3 = "workup3"
        lcms = "lcms"
        xchem = "xchem"
        nmr = "nmr"
        startingmaterial = "startingmaterial"
        solvent = "solvent"

    actionsession_id = models.ForeignKey(ActionSession, on_delete=models.CASCADE)
    reaction_id = models.ForeignKey(
        Reaction, related_name="addactions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    fromplatetype = models.CharField(choices=PlateType.choices, max_length=20)
    toplatetype = models.CharField(choices=PlateType.choices, max_length=20)
    smiles = models.CharField(max_length=255)
    calcunit = models.CharField(
        choices=CalcUnit.choices, default="moleq", max_length=10
    )
    volume = models.FloatField()
    volumeunit = models.CharField(
        choices=VolumeUnit.choices, default="ul", max_length=2
    )
    molecularweight = models.FloatField()
    solvent = models.CharField(max_length=255, null=True)
    concentration = models.FloatField(null=True)


class ExtractAction(models.Model):
    """Django model to define an extract action

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    fromplatetype: CharField
        The plate the add action is moving material from
    toplatetype: CharField
        The plate the add action is moviong material to
    layer: CharField
        The layer to extract
    volume: FloatField
        The volume to extract
    volumeunit: CharField
        The unit of the volume being extracted (default=ul)
    molecularweight: FloatField
        The molecular weight of the compound being added
    bottomlayervolume:
        The volume of the bottom layer (ul)
    bottomlayervolumeunit:
        The otional unit used to calculate the volume of the bottom layer
    solvent: CharField
        Optional solvent used to dilute the material being added
    concentration: FloatField
        Optional concentration of the material solution prepared
    """

    class Layer(models.TextChoices):
        top = "top"
        bottom = "bottom"

    class PlateType(models.TextChoices):
        reaction = "reaction"
        workup1 = "workup1"
        workup2 = "workup2"
        workup3 = "workup3"
        lcms = "lcms"
        xchem = "xchem"
        nmr = "nmr"
        startingmaterial = "startingmaterial"
        solvent = "solvent"

    actionsession_id = models.ForeignKey(ActionSession, on_delete=models.CASCADE)
    reaction_id = models.ForeignKey(
        Reaction, related_name="extractactions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    fromplatetype = models.CharField(choices=PlateType.choices, max_length=20)
    toplatetype = models.CharField(choices=PlateType.choices, max_length=20)
    layer = models.CharField(choices=Layer.choices, default="bottom", max_length=10)
    smiles = models.CharField(max_length=255)
    volume = models.FloatField()
    volumeunit = models.CharField(default="ul", max_length=2)
    molecularweight = models.FloatField()
    bottomlayervolume = models.FloatField(null=True)
    bottomlayervolumeunit = models.CharField(default="ul", max_length=2)
    solvent = models.CharField(max_length=255, null=True)
    concentration = models.FloatField(null=True)


class MixAction(models.Model):
    """Django model to define a mix action

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    platetype: CharField
        The plate being mixed
    repetitions: IntField
        The number of mixes
    volume: FloatField
        The volume to mix
    volumeunit: CharField
        The unit of the volume being mixed (default=ul)
    """

    class CalcUnit(models.TextChoices):
        masseq = "masseq"
        ul = "ul"

    class VolumeUnit(models.TextChoices):
        ul = "ul"
        ml = "ml"

    class PlateType(models.TextChoices):
        reaction = "reaction"
        workup1 = "workup1"
        workup2 = "workup2"
        workup3 = "workup3"
        lcms = "lcms"
        xchem = "xchem"
        nmr = "nmr"
        startingmaterial = "startingmaterial"
        solvent = "solvent"

    actionsession_id = models.ForeignKey(ActionSession, on_delete=models.CASCADE)
    reaction_id = models.ForeignKey(
        Reaction, related_name="mixactions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    platetype = models.CharField(choices=PlateType.choices, max_length=20)
    repetitions = models.IntegerField()
    calcunit = models.CharField(
        choices=CalcUnit.choices, default="moleq", max_length=10
    )
    volume = models.FloatField()
    volumeunit = models.CharField(
        choices=VolumeUnit.choices, default="ul", max_length=2
    )


class StirAction(models.Model):
    """Django model to define a StirAction - the stir action details

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    number: IntegerField
        The number of the action to be executed in a list of action numbers
    duration: FloatField
        The duration of the stir action
    durationunit: CharField
        The duration unit of the stir action (default=hours)
    temperature: IntegerField
        The temperature of the stir action (default=25)
    temperatureunit: CharField
        The temperature unit of the stir action (default=degC)
    stirringspeed: CharField
        The speed of the stir action (default=normal)
    """

    class PlateType(models.TextChoices):
        reaction = "reaction"
        workup1 = "workup1"
        workup2 = "workup2"
        workup3 = "workup3"
        lcms = "lcms"
        xchem = "xchem"
        nmr = "nmr"
        startingmaterial = "startingmaterial"
        solvent = "solvent"

    class TemperatureUnit(models.TextChoices):
        degcel = "degC"
        kelvin = "K"

    class Unit(models.TextChoices):
        seconds = "seconds"
        minutes = "minutes"
        hours = "hours"

    class Speed(models.TextChoices):
        gentle = "gentle"
        normal = "normal"
        vigorous = "vigorous"

    actionsession_id = models.ForeignKey(
        ActionSession, related_name="stiractions", on_delete=models.CASCADE
    )
    reaction_id = models.ForeignKey(
        Reaction, related_name="stiractions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    platetype = models.CharField(choices=PlateType.choices, max_length=20)
    duration = models.FloatField()
    durationunit = models.CharField(
        choices=Unit.choices, default=Unit.hours, max_length=10
    )
    temperature = models.IntegerField()
    temperatureunit = models.CharField(
        choices=TemperatureUnit.choices, default=TemperatureUnit.degcel, max_length=10
    )
    stirringspeed = models.CharField(
        choices=Speed.choices, default=Speed.normal, max_length=10
    )


# Models for capturing OT session, Deck, Plates and Wells
class OTProject(models.Model):
    """Django model to define an OTProject - an OT project will
       have one or more batch protocols for a project

    Parameters
    ----------
    project_id: ForeignKey
        Foreign key linking an OT protocol action to a reaction
    init_date: DateTimeField
        The date the OT project was created (autofield)
    name: CharField
        The name of the OT project
    """

    project_id = models.ForeignKey(Project, on_delete=models.CASCADE)
    init_date = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=150)


class OTBatchProtocol(models.Model):
    """Django model to define an OTBatchProtocol - OT protocols for a batch of
       targets

    Parameters
    ----------
    otproject_id: ForeignKey
        Foreign key linking a batch OT protocol with a OT project
    batch_id: ForeignKey
        Foreign key linking a batch OT protocol with batch of targets
    celery_taskid: CharField
        The Celery task id when a new OT batch protocol is created
    zipfile: FileField
        File link to a zip file of all the OT protocols required for executing the
        synthesis of a batch of targets
    """

    otproject_id = models.ForeignKey(OTProject, on_delete=models.CASCADE)
    batch_id = models.ForeignKey(
        Batch, related_name="otbatchprotocols", on_delete=models.CASCADE
    )
    celery_taskid = models.CharField(max_length=50)
    zipfile = models.FileField(upload_to="otbatchprotocols/", max_length=255, null=True)


class OTSession(models.Model):
    """Django model to define an OT Session - a OT session is a session (Reaction, Analysis)
       that needs to be executed on the OT

    Parameters
    ----------
    otbatchprotocol_id: ForeignKey
        Foreign key linking an OT session to a OT batch protocol
    reactionstep: IntegerField
        The reaction step that the session is being executed for
    sessiontype: CharField
        The type of session ebing excecuted
    """

    class SessionType(models.TextChoices):
        reaction = "reaction"
        workup = "workup"
        lcmsprep = "lcmsprep"

    otbatchprotocol_id = models.ForeignKey(
        OTBatchProtocol, related_name="otsessions", on_delete=models.CASCADE
    )
    reactionstep = models.IntegerField()
    sessiontype = models.CharField(
        choices=SessionType.choices, default=SessionType.reaction, max_length=10
    )


class Deck(models.Model):
    """Django model to define a Deck Session - an OT Deck

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a deck to an OT session
    numberslots: IntegerField
        The number of deck slots (default=11)
    slotavailable: BooleanField
        If a slot is still available on the deck (default=True)
    indexslotavailable: IntegerField
        The index (1-11) of the deck slot available
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="otdecks", on_delete=models.CASCADE
    )
    numberslots = models.IntegerField(default=11)
    slotavailable = models.BooleanField(default=True)
    indexslotavailable = models.IntegerField(default=1)


class Pipette(models.Model):
    """Django model to define a Pipette - an OT Pipette

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a deck to an OT session
    position: CharField
        The position (right or left) of the pipette for the OT session
    labware: CharField
        The name of the OT labware eg. p300_single
    type: CharField
        The type of pipette could be single or multi channel
    maxvolume: FloatField
        The maximum volume (ul) of the pipette
    name: CharField
        The name of the pipette - this is a combination of the OT labware name (p300_single)
        and the position of the pipette uses in the session (right) eg. right_p300_single
    """

    class Position(models.TextChoices):
        right = "Right"
        left = "Left"

    otsession_id = models.ForeignKey(
        OTSession, related_name="pipettes", on_delete=models.CASCADE
    )
    position = models.CharField(choices=Position.choices, max_length=10)
    labware = models.CharField(max_length=100)
    type = models.CharField(max_length=255)
    maxvolume = models.FloatField(default=300)
    name = models.CharField(max_length=255)


class TipRack(models.Model):
    """Django model to define a TipRack - an OT tiprack

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a tiprack to an OT session
    deck_id: ForeignKey
        Foreign key linking a tiprack to an OT deck
    labware: CharField
        The name of the OT labware eg. opentrons_96_tiprack_300ul
    index: IntegerField
        The deck index (1-11) of the tiprack
    name: CharField
        The name of the tiprack - this is a combination of the OT labware name (opentrons_96_tiprack_300ul)
        and the deck index of the tiprack eg. opentrons_96_tiprack_300ul_2
    """

    otsession_id = models.ForeignKey(OTSession, on_delete=models.CASCADE)
    deck_id = models.ForeignKey(Deck, on_delete=models.CASCADE)
    labware = models.CharField(max_length=255)
    index = models.IntegerField()
    name = models.CharField(max_length=255)


class Plate(models.Model):
    """Django model to define a Plate - an OT plate

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    deck_id: ForeignKey
        Foreign key linking a plate to an OT deck
    labware: CharField
        The name of the OT labware eg. labcyte_384_wellplate_100ul
    index: IntegerField
        The deck index (1-11) of the plate
    name: CharField
        The name of the plate
    type: CharField
        The type of plate eg. analyse and reaction plate
    maxwellvolume: FloatField
        The maximum plate well volume (ul)
    numberwells: IntegerField
        The number of plate wells
    wellavailable: BooleanField
        If a well is available on a plate
    indexswellavailable: IntegerField
        The index of the well available. Wells are occupied in increasing
        index starting from indices: A1, B1, C1 or 0, 1, 2 etc
    """

    class PlateType(models.TextChoices):
        reaction = "reaction"
        workup1 = "workup1"
        workup2 = "workup2"
        workup3 = "workup3"
        lcms = "lcms"
        xchem = "xchem"
        nmr = "nmr"
        startingmaterial = "startingmaterial"
        solvent = "solvent"

    otsession_id = models.ForeignKey(
        OTSession, related_name="otplates", on_delete=models.CASCADE
    )
    deck_id = models.ForeignKey(Deck, on_delete=models.CASCADE)
    labware = models.CharField(max_length=255)
    index = models.IntegerField()
    name = models.CharField(max_length=255)
    type = models.CharField(choices=PlateType.choices, max_length=55, null=True)
    maxwellvolume = models.FloatField()
    numberwells = models.IntegerField()
    wellavailable = models.BooleanField(default=True)
    indexswellavailable = models.IntegerField(default=0)


class Well(models.Model):
    """Django model to define a Well - an OT plate we

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    plate_id: ForeignKey
        Foreign key linking a well to a plate
    method_id: ForeignKey
        Optional foreign key linking a well to a method
    reaction_id: ForeignKey
        Optional foreign key linking a well to a reaction
    index: IntegerField
        The deck index (1-11) of the plate
    type: CharField
        The type of well eg. analyse and reaction well
    volume: FloatField
        The optional volume of the contents in the well (ul)
    smiles: CharField
        The optional SMILES of the well contents
    concentration: FloatField
        The optional concentration of the well contents
    solvent: CharField
        The optional solvent used for diluting the well contents
    reactantfornextstep: BooleanField
        Wether the contents are used for the next reaction step (default=False)
    available: BooleanField
        If the well is available w.r.t still containing it's contents
        vs. being empty (default=True)
    """

    class WellType(models.TextChoices):
        reaction = "reaction"
        workup1 = "workup1"
        workup2 = "workup2"
        workup3 = "workup3"
        lcms = "lcms"
        xchem = "xchem"
        nmr = "nmr"
        startingmaterial = "startingmaterial"
        solvent = "solvent"

    otsession_id = models.ForeignKey(
        OTSession, related_name="otwells", on_delete=models.CASCADE
    )
    plate_id = models.ForeignKey(Plate, on_delete=models.CASCADE)
    method_id = models.ForeignKey(Method, on_delete=models.CASCADE, null=True)
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE, null=True)
    clonewellid = models.IntegerField(null=True)
    index = models.IntegerField()
    type = models.CharField(choices=WellType.choices, max_length=55)
    volume = models.FloatField(null=True)
    smiles = models.CharField(max_length=255, null=True)
    concentration = models.FloatField(null=True)
    solvent = models.CharField(max_length=255, null=True)
    reactantfornextstep = models.BooleanField(default=False)
    available = models.BooleanField(default=True)


class CompoundOrder(models.Model):
    """Django model to define a CompoundOrder - a csv
       file of ordering information that includes SMILES,
       plate name, well id, amount, solvent and concentration
       required for the syhtesis of the reaction step

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    ordercsv: FileField
        The csv file of the ordering information for executing
        a reaction step on the OpenTrons
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="compoundorders", on_delete=models.CASCADE
    )
    ordercsv = models.FileField(upload_to="compoundorders/", max_length=255)


class OTScript(models.Model):
    """Django model to define a OTScript - a Python script
       for executing a reaction OpenTrons protocol

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    otscript: FileField
        The Python OpenTrons protocol file
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="otscripts", on_delete=models.CASCADE
    )
    otscript = models.FileField(upload_to="otscripts/", max_length=255)


class SolventPrep(models.Model):
    """Django model to define a SolventPrep - a csv
       file detailing the prepreation required for diluting previous
       reaction step products for use in the next reaction

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    solventprepcsv: FileField
        The csv file with solvent amount (ul), plate name and well index
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="solventpreps", on_delete=models.CASCADE
    )
    solventprepcsv = models.FileField(upload_to="solventprep/", max_length=255)
