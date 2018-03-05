from django.test import TestCase
from loader.functions import desalt_compound,NeutraliseCharges,get_path_or_none,sanitize_mol

from loader.loaders import add_target,add_mol,add_comp,add_prot,analyse_mols,analyse_target

# TEST ALL THE ABOVE
class LoaderFunctionsTestCase(TestCase):

    def test_desalt(self):
        test_input = "CCCCCC.Cl"
        test_output = "CCCCCC"
        output = desalt_compound(test_input)
        self.assertEqual(output,test_output)


    def test_neutralise(self):
        test_input = "CCC[NH+]CCC.Cl"
        test_output = ("CCCNCCC.Cl",True)
        output = NeutraliseCharges(test_input)
        self.assertEqual(output,test_output)

    def test_sanitize_mol(self):
        test_input = "CCC[NH+]CCC.Cl"
        test_output = "CCCNCCC"
        output = sanitize_mol(test_input)
        self.assertEqual(output,test_output)
