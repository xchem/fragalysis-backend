"""Checks validation of file for uploading to CAR"""
import math
from typing import BinaryIO
import inspect
import pandas as pd
from rdkit import Chem

from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import (
    canonSmiles,
    getAddtionOrder,
    checkReactantSMARTS,
    combiChem,
)

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
            expected_no_columns = 4
            expected_column_names = [
                "target-SMILES",
                "target-names",
                "amount-required-mg",
                "batch-tag",
            ]
            self.checkNumberColumns(
                expected_no_columns=expected_no_columns,
                expected_column_names=expected_column_names,
            )
            if self.validated:
                self.checkColumnNames(expected_column_names=expected_column_names)
            if self.validated:
                self.target_smiles = [
                    canonSmiles(smi.strip()) for smi in self.df["target-SMILES"]
                ]
                self.df["target-SMILES"] = self.target_smiles
                self.checkSMILES(
                    df_rows_index=self.index_df_rows,
                    smiles=self.target_smiles,
                    smiles_type="target",
                )
                if self.validated:
                    self.checkIsNumber()
                if self.validated:
                    self.checkIsString()

    def validateCustomChem(self):
        max_no_steps = max(self.df["no-steps"])
        reaction_numbers = list(range(1, max_no_steps + 1))
        expected_no_columns = (max_no_steps * 4) + 4
        expected_reactant_1_column_names = [
            "reactant-1-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reactant_2_column_names = [
            "reactant-2-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reaction_name_column_names = [
            "reaction-name-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_product_column_names = [
            "reaction-product-smiles-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_column_names = (
            [
                "target-names",
                "no-steps",
                "amount-required-mg",
                "batch-tag",
            ]
            + expected_reactant_1_column_names
            + expected_reactant_2_column_names
            + expected_product_column_names
            + expected_reaction_name_column_names
        )
        self.checkNumberColumns(
            expected_no_columns=expected_no_columns,
            expected_column_names=expected_column_names,
        )
        if self.validated:
            self.checkColumnNames(expected_column_names=expected_column_names)
        if self.validated:
            self.target_names = self.df["target-names"]
            self.batchtags = self.df["batch-tag"]
            self.amounts = self.df["amount-required-mg"]
            self.nosteps = self.df["no-steps"]
            self.target_smiles = []
            self.products = []
            self.reactant_pair_smiles = []
            self.reaction_names = []
            for index, row in self.df.iterrows():
                no_reaction_steps = row["no-steps"]
                reaction_numbers_group = list(range(1, no_reaction_steps + 1))
                for reaction_number in reaction_numbers_group:
                    reactant_1_SMILES = [
                        reactant.strip()
                        for reactant in row["reactant-1-{}".format(reaction_number)]
                        if str(reactant) != "nan"
                    ]

                    reactant_2_SMILES = [
                        reactant.strip()
                        for reactant in row["reactant-2-{}".format(reaction_number)]
                        if str(reactant) != "nan"
                    ]

                    reactant_pair_smiles = [reactant_1_SMILES, reactant_2_SMILES]
                    self.checkSMILES(
                        df_rows_index=index,
                        smiles=reactant_pair_smiles,
                        smiles_type="reactant_pair",
                    )
                    if self.validated:
                        (
                            reactant_pair_smiles_ordered,
                            product_smiles,
                        ) = self.checkReaction(
                            reactant_pair_smiles=reactant_pair_smiles,
                            reaction_names=row[
                                "reaction-name-{}".format(reaction_number)
                            ],
                            product_smiles=row[
                                "reaction-product-smiles-{}".format(reaction_number)
                            ],
                        )
                    reactant_pair_smiles_ordered = [
                        (canonSmiles(smi[0]), canonSmiles(smi[1]))
                        for smi in reactant_pair_smiles_ordered
                    ]

                    row[
                        "reaction-reactant-pair-smiles-{}".format(reaction_number)
                    ] = reactant_pair_smiles_ordered
                    if reaction_number == no_reaction_steps:
                        self.target_smiles = self.target_smiles + product_smiles

                reaction_names = list(
                    zip(
                        *[
                            row["reaction-name-{}".format(reactionnumber)]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                reactant_pair_smiles = list(
                    zip(
                        *[
                            row[
                                "reaction-reactant-pair-smiles-{}".format(
                                    reactionnumber
                                )
                            ]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                products = list(
                    zip(
                        *[
                            row["reaction-product-smiles-{}".format(reactionnumber)]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                self.products = self.products + products
                self.reactant_pair_smiles = (
                    self.reactant_pair_smiles + reactant_pair_smiles
                )
                self.reaction_names = self.reaction_names + reaction_names

            if self.validated:
                self.df = pd.DataFrame()
                self.df["batch-tag"] = self.batchtags
                self.df["target-names"] = self.target_names
                self.df["target-SMILES"] = self.target_smiles
                self.df["amount-required-mg"] = self.amounts
                self.df["no-steps"] = self.nosteps
                self.df["reactant-pair-smiles"] = self.reactant_pair_smiles
                self.df["reaction-name"] = self.reaction_names
                self.df["product-smiles"] = self.products
                self.checkIsNumber()

    def validateCustomCombiChem(self):
        max_no_steps = int(max(self.df["no-steps"]))
        reaction_numbers = list(range(1, max_no_steps + 1))
        expected_no_columns = (max_no_steps * 2) + 5
        expected_reactant_1_column_names = [
            "reactant-1-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reaction_name_column_names = [
            "reaction-name-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_column_names = (
            [
                "combi-group",
                "no-steps",
                "amount-required-mg",
                "batch-tag",
                "reactant-2-1",
            ]
            + expected_reactant_1_column_names
            + expected_reaction_name_column_names
        )
        self.checkNumberColumns(
            expected_no_columns=expected_no_columns,
            expected_column_names=expected_column_names,
        )
        if self.validated:
            self.checkColumnNames(expected_column_names=expected_column_names)

        if self.validated:
            self.target_names = []
            self.target_smiles = []
            self.nosteps = []
            self.products = []
            self.reactant_pair_smiles = []
            self.reaction_names = []
            self.batchtags = []
            self.amounts = []
            combi_grouped = self.df.groupby(["combi-group"])
            for combi_group_name, combi_group in combi_grouped:
                combi_group_info = {}
                combi_group = combi_group.reset_index()
                max_no_steps_combi_group = int(max(combi_group["no-steps"]))
                reaction_numbers_group = list(range(1, max_no_steps_combi_group + 1))
                columns_count = combi_group.nunique(
                    axis="rows", dropna=True
                )  # NB "nan" values (empty row values) not counted
                number_reactant_1s = [
                    columns_count["reactant-1-{}".format(reaction_number)]
                    for reaction_number in reaction_numbers_group
                    if "reactant-1-{}".format(reaction_number) in columns_count
                ]
                number_reactant_2s = [
                    columns_count["reactant-2-{}".format(reaction_number)]
                    for reaction_number in reaction_numbers_group
                    if "reactant-2-{}".format(reaction_number) in columns_count
                ]
                no_targets = int(math.prod(number_reactant_1s + number_reactant_2s))
                target_names = [
                    "{}-{}".format(combi_group_name, i) for i in range(no_targets)
                ]
                batch_tags = [combi_group.at[0, "batch-tag"]] * no_targets
                amounts = [combi_group.at[0, "amount-required-mg"]] * no_targets
                no_steps = [combi_group.at[0, "no-steps"]] * no_targets
                self.target_names = self.target_names + target_names
                self.batchtags = self.batchtags + batch_tags
                self.amounts = self.amounts + amounts
                self.nosteps = self.nosteps + no_steps
                for reaction_number in reaction_numbers_group:
                    reaction_combi_group_info = {}
                    if reaction_number == 1:
                        reactant_1_SMILES = [
                            reactant.strip()
                            for reactant in combi_group[
                                "reactant-1-{}".format(reaction_number)
                            ]
                            if str(reactant) != "nan"
                        ]

                        reactant_2_SMILES = [
                            reactant.strip()
                            for reactant in combi_group[
                                "reactant-2-{}".format(reaction_number)
                            ]
                            if str(reactant) != "nan"
                        ]

                    if reaction_number > 1:
                        reactant_1_SMILES = [
                            reactant.strip()
                            for reactant in combi_group[
                                "reactant-1-{}".format(reaction_number)
                            ]
                            if str(reactant) != "nan"
                        ]
                        reactant_2_SMILES = product_smiles[:number_reactant_pair_smiles]
                    reactant_pair_smiles = combiChem(
                        reactant_1_SMILES=reactant_1_SMILES,
                        reactant_2_SMILES=reactant_2_SMILES,
                    )
                    number_reactant_pair_smiles = len(reactant_pair_smiles)
                    if number_reactant_pair_smiles != no_targets:
                        reactant_pair_smiles = reactant_pair_smiles * (
                            no_targets // len(reactant_pair_smiles)
                        )

                    reaction_names = [
                        combi_group.at[0, "reaction-name-{}".format(reaction_number)]
                    ] * no_targets
                    reactant_pair_smiles_ordered, product_smiles = self.checkReaction(
                        reactant_pair_smiles=reactant_pair_smiles,
                        reaction_names=reaction_names,
                    )
                    reactant_pair_smiles_ordered = [
                        (canonSmiles(smi[0]), canonSmiles(smi[1]))
                        for smi in reactant_pair_smiles_ordered
                    ]
                    if reaction_number == max_no_steps_combi_group:
                        self.target_smiles = self.target_smiles + product_smiles
                    reaction_combi_group_info[
                        "reaction-name-{}".format(reaction_number)
                    ] = reaction_names
                    reaction_combi_group_info[
                        "reaction-reactant-pair-smiles-{}".format(reaction_number)
                    ] = reactant_pair_smiles_ordered
                    reaction_combi_group_info[
                        "reaction-product-smiles-{}".format(reaction_number)
                    ] = product_smiles
                    combi_group_info.update(reaction_combi_group_info)

                products = list(
                    zip(
                        *[
                            combi_group_info[
                                "reaction-product-smiles-{}".format(reactionnumber)
                            ]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                reactant_pair_smiles = list(
                    zip(
                        *[
                            combi_group_info[
                                "reaction-reactant-pair-smiles-{}".format(
                                    reactionnumber
                                )
                            ]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                reaction_names = list(
                    zip(
                        *[
                            combi_group_info["reaction-name-{}".format(reactionnumber)]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                self.products = self.products + products
                self.reactant_pair_smiles = (
                    self.reactant_pair_smiles + reactant_pair_smiles
                )
                self.reaction_names = self.reaction_names + reaction_names

            if self.validated:
                self.df = pd.DataFrame()
                self.df["batch-tag"] = self.batchtags
                self.df["target-names"] = self.target_names
                self.df["target-SMILES"] = self.target_smiles
                self.df["amount-required-mg"] = self.amounts
                self.df["no-steps"] = self.nosteps
                self.df["reactant-pair-smiles"] = self.reactant_pair_smiles
                self.df["reaction-name"] = self.reaction_names
                self.df["product-smiles"] = self.products
                self.checkIsNumber()

    def add_warning(self, field, warning_string):
        self.validate_dict["field"].append(field)
        self.validate_dict["warning_string"].append(warning_string)

    def checkColumnNames(self, expected_column_names):
        try:
            if not set(self.df_columns) == set(expected_column_names):
                self.add_warning(
                    field="name_columns",
                    warning_string="Column names should be set to: {}".format(
                        expected_column_names
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

    def checkNumberColumns(self, expected_no_columns, expected_column_names):
        try:
            if self.no_df_columns != expected_no_columns:
                self.add_warning(
                    field="number_columns",
                    warning_string="Found {} columns. Expected {} columns. Set and name columns to {} only".format(
                        self.no_df_columns,
                        expected_no_columns,
                        expected_column_names,
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

    def checkSMILES(
        self, df_rows_index: list[int], smiles: list[str], smiles_type: str
    ):
        """Checks if input SMILES from df is valid

        Parameters
        ----------
        df_rows_index: list[int]
            The index of the df rows being tested - for error reporting
        smiles: list[str]
            The SMILES being tested eg. Target or reactant pair SMILES
        smiles_type: str
            The type of SMILES being tested eg. target or reactant_pair

        """
        try:
            for index, smi in zip(df_rows_index, smiles):
                if all(isinstance(item, tuple) for item in smi):
                    mol_test = [Chem.MolFromSmiles(smi) for smi in smi]
                else:
                    mol_test = [Chem.MolFromSmiles(smi)]
                if None in mol_test:
                    self.add_warning(
                        field="check_smiles",
                        warning_string="Input {} smiles: '{}' at index {} is not a valid smiles".format(
                            smiles_type,
                            smi,
                            index,
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_smiles",
                warning_string="Input {} smiles check failed with error: {}".format(
                    smiles_type, e
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

    def checkReaction(
        self,
        reactant_pair_smiles: list,
        reaction_names: list[str],
        product_smiles: str = None,
    ):
        try:
            product_created_smiles = []
            reactant_pair_smiles_ordered = []
            for reactant_pair, reaction_name in zip(
                reactant_pair_smiles, reaction_names
            ):

                smarts = encoded_recipes[reaction_name]["recipes"]["standard"][
                    "reactionSMARTS"
                ]
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
                    if product_smiles:
                        product_mol = Chem.MolFromSmiles(product_smiles)
                    else:
                        product_mol = product_mols[0]
                    product_smi = Chem.MolToSmiles(product_mol)
                    reactant_smis = getAddtionOrder(
                        product_smi=product_smi,
                        reactant_SMILES=reactant_pair,
                        reaction_SMARTS=smarts,
                    )
                    product_created_smiles.append(product_smi)
                    reactant_pair_smiles_ordered.append(reactant_smis)
            return reactant_pair_smiles_ordered, product_created_smiles

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_reaction",
                warning_string="Reaction check failed with error: {}".format(e),
            )
            self.validated = False
