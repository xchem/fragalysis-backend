# from unittest import TestCase
# from rdkit import Chem
# from rdkit.Chem import AllChem

# from car.recipebuilder.encodedrecipes import encoded_recipes

# from car.utils import (
#     calculateMols,
#     canonSmiles,
#     combiChem,
#     createSVGString,
#     createReactionSVGString,
#     getAddtionOrder,
#     checkReactantSMARTS,
#     getChemicalName,
#     getPubChemCAS,
#     getPubChemCompound,
# )

# from car.tests.testdata.indata.testutils import (
#     snar_reactant_smiles_one,
#     snar_reactant_smiles_two,
# )
# from car.tests.testdata.outdata.testutils import (
#     snar_combo_equal,
#     snar_combo_unequal,
#     svg_str,
#     svg_reaction_str,
# )


# class ChemistryFunctionsTestCase(TestCase):
#     def setUp(self) -> None:
#         self.smiles = "C1COC2=C(C3=C(C(=C21)CCN)OCC3)Br"
#         self.snar_reactant_smiles_one = snar_reactant_smiles_one
#         self.snar_reactant_smiles_two = snar_reactant_smiles_two
#         self.snar_reactant_smiles_tuple = (
#             snar_reactant_smiles_one[0],
#             snar_reactant_smiles_two[0],
#         )
#         self.snar_encoded_smarts = encoded_recipes[
#             "N-nucleophilic aromatic substitution"
#         ]["reactionSMARTS"]
#         self.snar_product_smiles = "O=C(O)Cc1ccc(Nc2ccccc2)cc1F"
#         self.snar_product_mols = Chem.MolFromSmiles(self.snar_product_smiles)
#         self.svg_str = self.strip_white_space(str=svg_str)
#         self.reaction_smarts = AllChem.ReactionFromSmarts(
#             "{}.{}>>{}".format(
#                 self.snar_reactant_smiles_one[0],
#                 self.snar_reactant_smiles_two[0],
#                 self.snar_product_smiles,
#             )
#         )
#         self.reaction_svg_str = self.strip_white_space(svg_reaction_str)

#     def strip_white_space(self, str):
#         return str.replace(" ", "").replace("\t", "").replace("\n", "")

#     def test_calculate_product_mols(self):
#         target_mass = 10
#         test_product_mols = calculateMols(
#             target_mass=target_mass, target_SMILES=self.smiles
#         )

#         self.assertAlmostEqual(
#             first=test_product_mols,
#             second=3.533309327614245e-05,
#             places=4,
#             msg="incorrect product mols calculated",
#         )

#     def test_canon_smiles_correct(self):
#         test_canon_smiles = canonSmiles(smiles=self.smiles)
#         self.assertEqual(
#             test_canon_smiles,
#             "NCCc1c2c(c(Br)c3c1OCC3)OCC2",
#             "incorrect canonicalisation of SMILES",
#         )

#     def test_canon_smiles_incorrect(self):
#         test_canon_smiles = canonSmiles(smiles="OT Chemistry is possible")
#         self.assertEqual(
#             test_canon_smiles,
#             False,
#             "incorrect capture of bad SMILES input",
#         )

#     def test_combi_chem_equal(self):
#         test_all_possible_combinations = combiChem(
#             reactant_1_SMILES=self.snar_reactant_smiles_one,
#             reactant_2_SMILES=self.snar_reactant_smiles_two,
#         )
#         self.assertEqual(
#             test_all_possible_combinations,
#             snar_combo_equal,
#             "incorrect combinatorial product of two equal length lists of smiles",
#         )

#     def test_combi_chem_unequal(self):
#         test_all_possible_combinations = combiChem(
#             reactant_1_SMILES=self.snar_reactant_smiles_one[0:1],
#             reactant_2_SMILES=self.snar_reactant_smiles_two,
#         )
#         self.assertEqual(
#             test_all_possible_combinations,
#             snar_combo_unequal,
#             "incorrect combinatorial product of two non-equal length lists of smiles",
#         )

#     def test_create_svg_string(self):
#         test_svg_str = createSVGString(smiles=self.smiles)
#         test_svg_str = self.strip_white_space(test_svg_str)
#         self.assertEqual(
#             test_svg_str,
#             self.svg_str,
#             "incorrect creation of a svg string from SMILES",
#         )

#     def test_create_reaction_svg_string(self):
#         test_svg_str = createReactionSVGString(smarts=self.reaction_smarts)
#         test_svg_str = self.strip_white_space(test_svg_str)
#         self.assertEqual(
#             test_svg_str,
#             self.reaction_svg_str,
#             "incorrect creation of a reaction svg string from reaction SMARTS",
#         )

#     def test_get_addition_order_success(self):
#         test_ordered_smis = getAddtionOrder(
#             product_smi=self.snar_product_smiles,
#             reactant_SMILES=self.snar_reactant_smiles_tuple,
#             reaction_SMARTS=self.snar_encoded_smarts,
#         )
#         self.assertEqual(
#             test_ordered_smis,
#             list(self.snar_reactant_smiles_tuple),
#             "incorrect addtion order for a SNAr encoded recipe SMARTS",
#         )

#     def test_get_addition_order_fail(self):
#         first_reactant_smiles = "OT Chemistry is possible"
#         second_reactant_smiles = self.snar_reactant_smiles_one[0]
#         reactant_SMILES = (first_reactant_smiles, second_reactant_smiles)
#         test_ordered_smis = getAddtionOrder(
#             product_smi=self.snar_product_smiles,
#             reactant_SMILES=reactant_SMILES,
#             reaction_SMARTS=self.snar_encoded_smarts,
#         )
#         self.assertEqual(
#             test_ordered_smis,
#             None,
#             "incorrect SMILES input should return None for addtion order",
#         )

#     def test_check_reactant_smarts_success(self):
#         test_product_mols = checkReactantSMARTS(
#             reactant_SMILES=self.snar_reactant_smiles_tuple,
#             reaction_SMARTS=self.snar_encoded_smarts,
#         )

#         self.assertEqual(
#             len(test_product_mols),
#             2,
#             "incorrect length of product mols for testing reaction SMARTS",
#         )
#         self.assertEqual(
#             test_product_mols[1],
#             self.snar_product_smiles,
#             "incorrect product SMILES match for testing reaction SMARTS",
#         )


# class PubChemFunctionsTestCase(TestCase):
#     def setUp(self) -> None:
#         self.smiles = "C1COC2=C(C3=C(C(=C21)CCN)OCC3)Br"
#         self.compound = getPubChemCompound(smiles=self.smiles)

#     def test_get_pubchem_compound_success(self):
#         test_compound = getPubChemCompound(smiles=self.smiles)
#         self.assertEqual(
#             test_compound.cid,
#             10265873,
#             "incorrect PubChem id for getting a compound from PubChem",
#         )
#         self.assertEqual(
#             test_compound.cas,
#             "733720-95-1",
#             "incorrect CAS number for PubChem compound search",
#         )

#     def test_get_pubchem_compound_fail(self):
#         test_compound = getPubChemCompound(smiles="OT chemistry is possible")
#         self.assertEqual(
#             test_compound,
#             None,
#             "PubChem search should yield None response",
#         )

#     def test_get_chemical_name_success(self):
#         test_name = getChemicalName(smiles=self.smiles)
#         self.assertEqual(
#             test_name,
#             "2-(4-bromo-2,3,6,7-tetrahydrofuro[2,3-f][1]benzofuran-8-yl)ethanamine",
#             "incorrect PubChem IUPAC name for compound",
#         )

#     def test_get_chemical_name_fail(self):
#         test_name = getChemicalName(smiles="OT chemistry is possible")
#         self.assertEqual(
#             test_name,
#             None,
#             "PubChem name search should fail and return None",
#         )

#     def test_get_pubchem_cas(self):
#         test_cas = getPubChemCAS(compound=self.compound)
#         self.assertEqual(test_cas, "733720-95-1", "incorrect CAS number returned")
