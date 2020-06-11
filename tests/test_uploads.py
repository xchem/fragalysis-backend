import unittest
from viewer.sdf_check import *

from rdkit import Chem
from rdkit.Chem import Descriptors


class ValidateTestValid(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sdf_file = '/code/tests/test_data/compund-set_AlessandroContini-NwatMM-GBSAcovalentrun1.sdf'
        cls.Target = 'Mpro'
        cls.version = 'ver_1.2'
        cls.suppl = Chem.SDMolSupplier(cls.sdf_file)
        cls.blank_mol = cls.suppl[0]
        cls.other_mols = cls.suppl[1:]
        cls.validate_dict = {'molecule_name': [],
                             'field': [],
                             'warning_string': []}

    @classmethod
    def tearDownClass(cls):
        pass

    def test_check_mol_props(self):
        expected_dict = {'molecule_name': [],
                         'field': [],
                         'warning_string': []}
        for mol in self.other_mols:
            check_mol_props(mol, self.validate_dict)

        self.assertEqual(self.validate_dict, expected_dict)

    def test_check_missing_field(self):
        pass

    def test_check_name_characters(self):
        pass

    def test_check_url(self):
        pass

    def test_check_field_populated(self):
        pass

    def test_check_blank_prop(self):
        pass

    def test_check_blank_mol_props(self):
        pass

    def test_check_ver_name(self):
        pass

    def test_check_SMILES(self):
        pass

    def test_check_pdb(self):
        pass

    def test_check_refmol(self):
        pass

    def test_check_sdf(self):
        pass

    def test_add_warning(self):
        pass

    def test_check_compound_set(self):
        pass

    def test_check_property_descriptions(self):
        pass
