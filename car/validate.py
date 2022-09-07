"""Checks validation of file for uploading to CAR"""
from typing import BinaryIO
import inspect
import pandas as pd
from rdkit import Chem

from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import canonSmiles, getAddtionOrder, checkReactantSMARTS, combiChem

import logging

logger = logging.getLogger(__name__)


class ValidateFile(object):
    """
    Creates a validate object for checking file validation for upload
    """

    def __init__(self, csv_to_validate: BinaryIO, validate_type: str):
        """ValidateFile constructor

        Parameters
        ----------

        csv_to_validate: IO
            The uploaded csv file for validating
        """
        self.df = pd.read_csv(csv_to_validate, encoding="utf8")
        self.df_columns = self.df.columns
        self.no_df_columns = len(self.df_columns)
        self.index_df_rows = range(0, len(self.df), 1)
        self.upload_type = validate_type
        self.validate_dict = {"field": [], "warning_string": []}
        self.validated = True

        if self.upload_type == "custom-chem":
            self.validateCustomChem()
        if self.upload_type == "combi-custom-chem":
            self.validateCustomCombiChem()

        if self.upload_type == "retro-API":
            self.expected_no_columns = 3
            self.expected_column_names = ["targets", "amount-required-mg", "batch-tag"]
            self.checkNumberColumns()
            if self.validated:
                self.checkColumnNames()
            if self.validated:
                self.target_smiles = [
                    canonSmiles(smi.strip()) for smi in self.df["targets"]
                ]
                self.df["targets"] = self.target_smiles
                self.checkTargetSMILES()
                if self.validated:
                    self.checkIsNumber()
                if self.validated:
                    self.checkIsString()

    def validateCustomChem(self):
        """Validates the csv file if the custom chemistry
        option selected
        """
        self.expected_no_columns = 5
        self.expected_column_names = [
            "reactant-1",
            "reactant-2",
            "reaction-name",
            "amount-required-mg",
            "batch-tag",
        ]
        self.checkNumberColumns()
        if self.validated:
            self.checkColumnNames()
        if self.validated:
            self.reactant_pair_smiles = [
                reactants
                for reactants in zip(self.df["reactant-1"], self.df["reactant-2"])
            ]
            self.df["reactant-pair-smiles"] = self.reactant_pair_smiles
            self.checkReactantSMILES()
            if self.validated:
                self.reactant_pair_smiles = [
                    (canonSmiles(smi[0]), canonSmiles(smi[1]))
                    for smi in self.reactant_pair_smiles
                ]
                self.reaction_names = self.df["reaction-name"]
                self.checkReaction()
                if self.validated:
                    self.df["reactant-pair-smiles"] = self.reactant_pair_smiles_ordered
                    self.df["target-smiles"] = self.product_smiles
                    self.checkIsNumber()
                if self.validated:
                    self.checkIsString()

    def validateCustomCombiChem(self):
        max_no_steps = max(self.df["no-steps"])
        self.expected_no_columns = (max_no_steps * 3) + 2
        reactant_1_column_names = []  # Must finish column naming!
        reactant_2_column_names = []

        self.expected_column_names = [
            "reactant-1",
            "reactant-2",
            "reaction-name",
            "amount-required-mg",
            "batch-tag",
        ]
        self.checkNumberColumns()
        # if self.validated:
        #     self.checkColumnNames()

        if self.validated:
            self.reactant_pair_smiles = []
            self.reaction_names = []
            self.batchtags = []

            combi_grouped = self.df.groupby(["combi-group"])
            for combi_group in combi_grouped:
                max_no_steps_combi_group = max(combi_group["no-steps"])
                for i in range(max_no_steps_combi_group):
                    reaction_step_df = combi_group.filter(regex="-{}$".format(i))
                    reaction_tag_grouped = reaction_step_df.groupby(
                        ["reaction-name-{}".format(i), "batch-tag"]
                    )
                    for name, group in reaction_tag_grouped:
                        group = group.reset_index()
                        reactant_1_SMILES = set(
                            [
                                reactant
                                for reactant in group["reactant-1-{}".format(i)]
                                if str(reactant) != "nan"
                            ]
                        )
                        reactant_2_SMILES = set(
                            [
                                reactant
                                for reactant in group["reactant-2-{}".format(i)]
                                if str(reactant) != "nan"
                            ]
                        )
                        reactant_pair_smiles = combiChem(
                            reactant_1_SMILES=reactant_1_SMILES,
                            reactant_2_SMILES=reactant_2_SMILES,
                        )
                        reaction_names = [name[0]] * len(reactant_pair_smiles)
                        batchtags = [group.at[0, "batch-tag"]] * len(
                            reactant_pair_smiles
                        )
                        self.reactant_pair_smiles = (
                            self.reactant_pair_smiles + reactant_pair_smiles
                        )
                        self.reaction_names = self.reaction_names + reaction_names
                        self.batchtags = self.batchtags + batchtags

                    self.checkReactantSMILES()
                    if self.validated:
                        self.reactant_pair_smiles = [
                            (canonSmiles(smi[0]), canonSmiles(smi[1]))
                            for smi in self.reactant_pair_smiles
                        ]
                        self.checkReaction()
                        if self.validated:
                            amount_required_mg = self.df.at[0, "amount-required-mg"]
                            self.df = pd.DataFrame()
                            self.df[
                                "reactant-pair-smiles"
                            ] = self.reactant_pair_smiles_ordered
                            self.df["target-smiles"] = self.product_smiles
                            self.df["reaction-name"] = self.reaction_names
                            self.df["amount-required-mg"] = [amount_required_mg] * len(
                                self.reactant_pair_smiles
                            )
                            self.df["batch-tag"] = self.batchtags
                            self.checkIsNumber()

    def add_warning(self, field, warning_string):
        self.validate_dict["field"].append(field)
        self.validate_dict["warning_string"].append(warning_string)

    def checkColumnNames(self):
        try:
            if not all(self.df_columns == self.expected_column_names):
                self.add_warning(
                    field="name_columns",
                    warning_string="Column names should be set to: {}".format(
                        self.expected_column_names
                    ),
                )
                self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="name_columns",
                warning_string="Column names check failed with error: {}".format(e),
            )
            self.validated = False

    def checkNumberColumns(self):
        try:
            if self.no_df_columns != self.expected_no_columns:
                self.add_warning(
                    field="number_columns",
                    warning_string="Found {} columns. Expected {} columns. Set and name columns to {} only".format(
                        self.no_df_columns,
                        self.expected_no_columns,
                        self.expected_column_names,
                    ),
                )
                self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="number_columns",
                warning_string="Number columns check failed with error: {}".format(e),
            )
            self.validated = False

    def checkTargetSMILES(self):
        try:
            for index, smi in zip(self.index_df_rows, self.target_smiles):
                canonsmiles = canonSmiles(smi)
                if not canonsmiles:
                    self.add_warning(
                        field="check_smiles",
                        warning_string="Input target smiles: '{}' at index {} is not a valid smiles".format(
                            smi,
                            index,
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_smiles",
                warning_string="Input target smiles check failed with error: {}".format(
                    e
                ),
            )
            self.validated = False

    def checkReactantSMILES(self):
        try:
            for index, smi_pair in zip(self.index_df_rows, self.reactant_pair_smiles):
                mols = [Chem.MolFromSmiles(smi) for smi in smi_pair]
                if None in mols:
                    none_test_indices = [
                        index for index, mol in enumerate(mols) if mol is None
                    ]
                    invalid_smiles = [smi_pair[index] for index in none_test_indices]
                    self.add_warning(
                        field="check_smiles",
                        warning_string="Input reactant smiles: ".join(
                            "{} ".format(*smi) for smi in invalid_smiles
                        )
                        + "at index {} is not a valid smiles".format(
                            index,
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_smiles",
                warning_string="Input reactant smiles check failed with error: {}".format(
                    e
                ),
            )
            self.validated = False

    def checkIsNumber(self):
        try:
            self.target_amounts = [amount for amount in self.df["amount-required-mg"]]
            for index, amount in zip(self.index_df_rows, self.target_amounts):
                if not isinstance(amount, (int, float)):
                    self.add_warning(
                        field="check_number",
                        warning_string="Target mass {} at index {} is not a valid number".format(
                            amount, index
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_number",
                warning_string="Input amount check failed with error: {}".format(e),
            )
            self.validated = False

    def checkIsString(self):
        try:
            self.batchtags = [tag.strip() for tag in self.df["batch-tag"]]
            for index, tag in zip(self.index_df_rows, self.batchtags):
                if not type(tag) == str:
                    self.add_warning(
                        field="check_string",
                        warning_string="Batch tag {} at index {} is not a valid string format".format(
                            tag, index
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_string",
                warning_string="Input batch tag check failed with error: {}".format(e),
            )
            self.validated = False

    def checkReaction(self):
        try:
            self.product_smiles = []
            self.reactant_pair_smiles_ordered = []
            no_reaction_tests = len(self.reaction_names)

            for index, reactant_pair, reaction_name in zip(
                range(no_reaction_tests), self.reactant_pair_smiles, self.reaction_names
            ):

                smarts = encoded_recipes[reaction_name]["reactionSMARTS"]
                product_mols = checkReactantSMARTS(
                    reactant_SMILES=reactant_pair, reaction_SMARTS=smarts
                )

                if not product_mols:
                    self.add_warning(
                        field="check_reaction",
                        warning_string="Reaction for reactants: {} and reaction: {} is not a valid reaction".format(
                            reactant_pair, reaction_name
                        ),
                    )
                    self.validated = False

                if product_mols:
                    product_mol = product_mols[0]
                    product_smi = Chem.MolToSmiles(product_mol)
                    reactant_smis = getAddtionOrder(
                        product_smi=product_smi,
                        reactant_SMILES=reactant_pair,
                        reaction_SMARTS=smarts,
                    )
                    self.product_smiles.append(product_smi)
                    self.reactant_pair_smiles_ordered.append(reactant_smis)
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_reaction",
                warning_string="Reaction check failed with error: {}".format(e),
            )
            self.validated = False
