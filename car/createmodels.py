"""Create Django models from IBM API output"""
from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from statistics import mean

from .mcule.apicalls import MCuleAPI

# Import standard models
from .models import (
    Project,
    Batch,
    Target,
    Method,
    Reaction,
    Product,
    Reactant,
    CatalogEntry,
    MculeQuote,
)

# Import action models
from .models import (
    AddAction,
    AnalyseAction,
    ExtractAction,
    FilterAction,
    QuenchAction,
    SetTemperatureAction,
    StirAction,
)

from .utils import calculateproductmols, createSVGString, createReactionSVGString


def createProjectModel(project_info):
    project = Project()
    project.name = project_info["projectname"]
    project.submittername = project_info["submittername"]
    project.submitterorganisation = project_info["submitterorganisation"]
    project.proteintarget = project_info["proteintarget"]
    project.save()
    return project.id, project.name


def createBatchModel(project_id, batch_tag, batch_id=None):
    batch = Batch()
    project_obj = Project.objects.get(id=project_id)
    batch.project_id = project_obj
    if batch_id:
        fk_batch_obj = Batch.objects.get(pk=batch_id)
        batch.batch_id = fk_batch_obj
    batch.batch_tag = batch_tag
    batch.save()
    return batch.id


def createTargetModel(batch_id, smiles, target_mass):
    """
    Function that creates a Target object
    if the csv file uploaded is validated and
    the user wants to upload the data

    project_id: string id of project created for upload
    smiles: string a valid smiles
    """
    target = Target()
    batch_obj = Batch.objects.get(id=batch_id)
    batch_tag = batch_obj.batch_tag
    target.batch_id = batch_obj
    target.smiles = smiles
    target.targetmols = calculateproductmols(target_mass, smiles)
    target.name = batch_tag
    target_svg_string = createSVGString(smiles)
    target_svg_fn = default_storage.save(
        "targetimages/" + target.name + ".svg", ContentFile(target_svg_string)
    )
    target.image = target_svg_fn
    target.targetmass = target_mass
    target.unit = "mg"
    target.save()
    return target.id


def createMethodModel(target_id, nosteps, otchem):
    method = Method()
    target_obj = Target.objects.get(id=target_id)
    method.target_id = target_obj
    method.nosteps = nosteps
    method.otchem = otchem
    method.save()

    return method.id


def createReactionModel(
    method_id, reaction_class, reaction_smarts, reaction_temperature=None
):
    reaction = Reaction()
    method_obj = Method.objects.get(id=method_id)
    reaction.method_id = method_obj
    reaction.reactionclass = reaction_class
    if reaction_temperature:
        reaction.reactiontemperature = reaction_temperature
    reaction_svg_string = createReactionSVGString(reaction_smarts)
    reaction_svg_fn = default_storage.save(
        "reactionimages/" + reaction_class + ".svg", ContentFile(reaction_svg_string)
    )
    reaction.reactionimage = reaction_svg_fn
    reaction.save()

    return reaction.id


def createProductModel(reaction_id, project_name, batch_tag, product_smiles):
    product = Product()
    product.name = "{}-{}".format(project_name, batch_tag)
    reaction_obj = Reaction.objects.get(id=reaction_id)
    product.reaction_id = reaction_obj
    product.smiles = product_smiles
    product_svg_string = createSVGString(product_smiles)
    product_svg_fn = default_storage.save(
        "productimages/" + product.name + ".svg", ContentFile(product_svg_string)
    )
    product.image = product_svg_fn
    product.save()


def createReactantModel(reaction_id, reactant_smiles):
    reactant = Reactant()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    reactant.reaction_id = reaction_obj
    reactant.smiles = reactant_smiles
    reactant.save()
    return reactant.id


def createCatalogEntryModel(catalog_entry, target_id=None, reactant_id=None):
    catalogentry = CatalogEntry()
    if target_id:
        target_obj = Target.objects.get(id=target_id)
        catalogentry.target_id = target_obj
    if reactant_id:
        reactant_obj = Reactant.objects.get(id=reactant_id)
        catalogentry.reactant_id = reactant_obj

    catalogentry.vendor = catalog_entry["catalogName"]
    catalogentry.catalogid = catalog_entry["catalogId"]

    if catalog_entry["catalogName"] == "generic":
        catalogentry.price = 0
        catalogentry.leadtime = 0

    if catalog_entry["purchaseInfo"]["isScreening"]:
        if catalog_entry["purchaseInfo"]["scrLeadTimeWeeks"] != "unknown":
            catalogentry.leadtime = catalog_entry["purchaseInfo"]["scrLeadTimeWeeks"]
        else:
            catalogentry.leadtime = None
        if catalog_entry["purchaseInfo"]["scrPriceRange"] != "unknown":
            priceinfo = catalog_entry["purchaseInfo"]["scrPriceRange"]
            catalogentry.priceinfo = priceinfo
            priceinfo = priceinfo.replace(" ", "")
            # Check type of range is less than or given range
            if priceinfo[0] == "<" or priceinfo[0] == ">":
                if "k" in priceinfo:
                    upperprice = int("".join(filter(str.isdigit, priceinfo))) * 1000
                else:
                    upperprice = int("".join(filter(str.isdigit, priceinfo)))
            if priceinfo[0] == "$":
                if "k" in priceinfo:
                    upperprice = (
                        int("".join(filter(str.isdigit, priceinfo.split("-")[1])))
                        * 1000
                    )
                else:
                    upperprice = int(
                        "".join(filter(str.isdigit, priceinfo.split("-")[1]))
                    )
            catalogentry.upperprice = upperprice
        else:
            catalogentry.price = None

    if not catalog_entry["purchaseInfo"]["isScreening"]:
        if catalog_entry["purchaseInfo"]["bbLeadTimeWeeks"] != "unknown":
            catalogentry.leadtime = catalog_entry["purchaseInfo"]["bbLeadTimeWeeks"]
        else:
            catalogentry.leadtime = None
        if catalog_entry["purchaseInfo"]["bbPriceRange"] != "unknown":
            priceinfo = catalog_entry["purchaseInfo"]["bbPriceRange"]
            catalogentry.priceinfo = priceinfo
            priceinfo = priceinfo.replace(" ", "")
            # Check type of range is less than or given range
            if priceinfo[0] == "<" or priceinfo[0] == ">":
                if "k" in priceinfo:
                    upperprice = int("".join(filter(str.isdigit, priceinfo))) * 1000
                else:
                    upperprice = int("".join(filter(str.isdigit, priceinfo)))
            if priceinfo[0] == "$":
                if "k" in priceinfo:
                    upperprice = (
                        int("".join(filter(str.isdigit, priceinfo.split("-")[1])))
                        * 1000
                    )
                else:
                    upperprice = int(
                        "".join(filter(str.isdigit, priceinfo.split("-")[1]))
                    )
            catalogentry.upperprice = upperprice
        else:
            catalogentry.price = None
    catalogentry.save()


class CreateEncodedActionModels(object):
    """
    Creates a createEncodedActionModels object for creating action models
    for a reaction
    """

    def __init__(
        self,
        actions: list,
        target_id: int,
        reaction_id: int,
        reactant_pair_smiles: list,
        reaction_name: str,
    ):
        """
        ValidateFile constructor
        Args:
            actions (list): List of actions
            project_id (int): Project model id
            reaction_id (int): Reaction model id for actions
            target_id (int): Target model id for reaction
            reactant_pair_smiles (list): List of reactant smiles
            reaction_name (str): Reaction name
        """
        self.mculeapi = MCuleAPI()
        self.actions = actions
        self.reaction_id = reaction_id
        self.reaction_obj = Reaction.objects.get(id=reaction_id)
        self.reactant_pair_smiles = reactant_pair_smiles
        self.reaction_name = reaction_name
        self.target_mols = Target.objects.get(id=target_id).targetmols
        self.mculeidlist = []
        self.amountslist = []

        for action in self.actions:
            self.createEncodedActionModel(action)

    def createEncodedActionModel(self, action):
        actionMethods = {
            "add": self.createAddActionModel,
            "analyse": self.createAnalyseActionModel,
            "extract": self.createExtractActionModel,
            "filter": self.createFilterActionModel,
            "quench": self.createQuenchActionModel,
            "set-temperature": self.createSetTemperatureActionModel,
            "stir": self.createStirActionModel,
        }

        action_type = action["name"]

        if action_type in actionMethods:
            actionMethods[action_type](action_type, action)
            return True
        else:
            logger.info(action_type)
            print(action)

    def calculateVolume(
        self, molar_eqv, conc_reagents=None, reactant_density=None, reactant_MW=None
    ):
        # NB need addition_order added to ncoded recipes - can we rather use SA score?
        mol_material = molar_eqv * self.target_mols

        if reactant_density:
            vol_material = ((mol_material * reactant_MW) / reactant_density) * 1e3
        else:
            vol_material = (mol_material / conc_reagents) * 1e6  # in uL
        return vol_material

    def calculateAmount(self, molar_eqv: float, reactant_MW: float):
        """ "
        Calculates amount of compound needed in mg

        Args:
            molar_eq (float): Molar equivalents required
            reactant_MW (float): Molecular weight of compound
        Returns:
            amount (float): Amount required in mg
        """
        mol_material = molar_eqv * self.target_mols
        mass_material = (mol_material / reactant_MW) * 1e6
        return mass_material

    def createAnalyseActionModel(self):
        analyse = AnalyseAction()
        analyse.reaction_id = self.reaction_obj
        analyse.actiontype = "analyse"
        analyse.actionno = self.action_no
        analyse.save()

    def createAddActionModel(self, action_type, action):
        try:
            if action["content"]["material"]["SMARTS"]:
                reactant_SMILES = self.reactant_pair_smiles[0]
                del self.reactant_pair_smiles[0]
                print(self.reactant_pair_smiles)
            if action["content"]["material"]["SMILES"]:
                reactant_SMILES = action["content"]["material"]["SMILES"]
            action_no = action["content"]["action_no"]
            molar_eqv = action["content"]["material"]["quantity"]["value"]
            concentration = action["content"]["material"]["concentration"]
            if not concentration:
                concentration = 0
            solvent = action["content"]["material"]["solvent"]
            mol = Chem.MolFromSmiles(reactant_SMILES)
            reactant_MW = Descriptors.ExactMolWt(mol)
            amount = self.calculateAmount(molar_eqv=molar_eqv, reactant_MW=reactant_MW)

            add = AddAction()
            add.reaction_id = self.reaction_obj
            add.actiontype = action_type
            add.actionno = action_no
            # material = getChemicalName(reactant_SMILES)
            material = reactant_SMILES
            add.material = material
            # if not material:
            #     add.material = str(self.reaction_id) + str(action_no) + "-" + reactant_SMILES
            # else:
            #     add.material = material
            mol = Chem.MolFromSmiles(reactant_SMILES)
            molecular_weight = Descriptors.ExactMolWt(mol)
            add.materialsmiles = reactant_SMILES
            # mculeinfo = self.mculeapi.getMCuleInfo(smiles=reactant_SMILES)
            # if mculeinfo:
            #     mculeid = mculeinfo[0]
            #     self.mculeidlist.append(mculeid)
            #     self.amountslist.append(amount)
            #     add.mculeid = mculeid
            #     add.mculeurl = mculeinfo[1]
            #     priceinfo = self.mculeapi.getMCulePrice(mculeid=mculeid, amount=amount)
            #     if priceinfo:
            #         add.mculeprice = priceinfo[0]
            #         add.mculedeliverytime = priceinfo[1]
            add.molecularweight = molecular_weight
            add_svg_string = createSVGString(reactant_SMILES)
            add_svg_fn = default_storage.save(
                "addactionimages/{}-{}-{}.svg".format(
                    self.reaction_id, action_no, material
                ),
                ContentFile(add_svg_string),
            )
            add.materialimage = add_svg_fn  # need material (common name)
            add.atmosphere = "air"

            if not solvent:
                reactant_density = action["content"]["material"]["density"]
                add.materialquantity = self.calculateVolume(
                    molar_eqv=molar_eqv,
                    reactant_density=reactant_density,
                    reactant_MW=reactant_MW,
                )

            if solvent:
                add.materialquantity = self.calculateVolume(
                    molar_eqv=molar_eqv,
                    conc_reagents=concentration,
                )
                add.solvent = solvent
            add.concentration = concentration
            add.save()

        except Exception as error:
            print(error)
            print(action)

    def createExtractActionModel(self, action_type, action):
        try:
            solvent = action["content"]["solvent"]["value"]
            quantity = action["content"]["solvent"]["quantity"]["value"]
            unit = action["content"]["solvent"]["quantity"]["unit"]
            repetitions = action["content"]["repetitions"]["value"]

            extract = ExtractAction()
            extract.reaction_id = self.reaction_obj
            extract.actiontype = action_type
            extract.actionno = self.action_no
            extract.solvent = solvent
            extract.solventquantity = quantity
            extract.solventquantityunit = unit
            extract.numberofrepetitions = repetitions
            extract.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createFilterActionModel(self, action_type, action):
        try:
            phasetokeep = action["content"]["phase_to_keep"]["value"]
            rinsingsolvent = action["content"]["rinsing_solvent"]["value"]
            rinsingsolventquantity = action["content"]["rinsing_solvent"]["quantity"][
                "value"
            ]
            rinsingsolventquantityunit = action["content"]["rinsing_solvent"][
                "quantity"
            ]["unit"]
            extractionforprecipitatesolvent = action["content"]["extraction_solvent"][
                "value"
            ]
            extractionforprecipitatesolventquantity = action["content"][
                "extraction_solvent"
            ]["quantity"]["value"]
            extractionforprecipitatesolventquantityunit = action["content"][
                "extraction_solvent"
            ]["quantity"]["unit"]

            filteraction = FilterAction()
            filteraction.reaction_id = self.reaction_obj
            filteraction.actiontype = action_type
            filteraction.actionno = self.action_no
            filteraction.phasetokeep = phasetokeep
            filteraction.rinsingsolvent = rinsingsolvent
            filteraction.rinsingsolventquantity = rinsingsolventquantity
            filteraction.rinsingsolventquantityunit = rinsingsolventquantityunit
            filteraction.extractionforprecipitatesolvent = (
                extractionforprecipitatesolvent
            )
            filteraction.extractionforprecipitatesolventquantity = (
                extractionforprecipitatesolventquantity
            )
            filteraction.extractionforprecipitatesolventquantityunit = (
                extractionforprecipitatesolventquantityunit
            )
            filteraction.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createQuenchActionModel(self, action_type, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            temperature = action["content"]["temperature"]

            quench = QuenchAction()
            quench.reaction_id = self.reaction_obj
            quench.actiontype = action_type
            quench.actionno = self.action_no
            quench.material = material
            quench.materialquantity = materialquantity
            quench.materialquantityunit = materialquantityunit
            if temperature:
                quench.temperature = temperature["value"]
            quench.dropwise = dropwise
            quench.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createSetTemperatureActionModel(self, action_type, action):
        try:
            temperature = action["content"]["temperature"]["value"]

            temp = SetTemperatureAction()
            temp.reaction_id = self.reaction_obj
            temp.actiontype = action_type
            temp.actionno = self.action_no
            temp.temperature = temperature
            temp.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createStirActionModel(self, action_type, action):
        try:
            action_no = action["content"]["action_no"]
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]["value"]
            temperatureunit = action["content"]["temperature"]["unit"]

            stir = StirAction()
            stir.reaction_id = self.reaction_obj
            stir.actiontype = action_type
            stir.actionno = action_no
            stir.duration = duration
            stir.durationunit = durationunit
            stir.temperature = temperature
            stir.temperatureunit = temperatureunit
            stir.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)


class CreateMculeQuoteModel(object):
    """
    Creates a CreateMculeQuoteModel object for creating a Mcule quote
    for a project
    """

    def __init__(
        self,
        mculeids: list,
        amounts: list,
        project_id: int,
    ):
        """
        ValidateFile constructor
        Args:
            mculeids (list): List of mcule ids
            project_id (int): Project model id
        """
        self.mculeidlist = [item for sublist in mculeids for item in sublist]
        self.amounts = amounts
        self.amountaverage = self.getAmountAverage()
        self.project_obj = Project.objects.get(id=project_id)
        self.mculeapi = MCuleAPI()
        self.createMculeQuoteModel()

    def getAmountAverage(self):
        if len(self.amounts) == 0:
            return 0
        if len(self.amounts) == 1:
            return self.amounts[0]
        else:
            mean([item for sublist in self.amounts for item in sublist])
            return mean

    def createMculeQuoteModel(self):
        quote_info = self.mculeapi.getTotalQuote(
            mculeids=self.mculeidlist, amount=self.amountaverage
        )

        if quote_info:
            try:
                quote = MculeQuote()
                quote.project_id = self.project_obj
                quote.quoteid = quote_info["quoteid"]
                quote.quoteurl = quote_info["quoteurl"]
                quote.quotecost = quote_info["quotecost"]
                quote.quotevaliduntil = quote_info["quotevaliduntil"]
                quote.save()

            except Exception as error:
                print(error)
