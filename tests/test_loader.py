from django.test import TestCase
from loader.functions import desalt_compound,NeutraliseCharges,sanitize_mol
from loader.loaders import add_target, add_mol, add_comp, add_prot, analyse_mols, analyse_target
from rdkit import Chem
from viewer.models import Molecule,Compound,Protein,Target
from scoring.models import MolGroup

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
        test_input = Chem.MolFromSmiles("CCC[NH+]CCC.Cl")
        test_output = "CCCNCCC"
        output = sanitize_mol(test_input)
        self.assertEqual(Chem.MolToSmiles(output),test_output)


class LoderLoaderTestCase(TestCase):

    def setUp(self):
        out_f = open('DUMMY_NEW.pdb', "w")
        out_f = open("DUMMY_NEW_TWO.pdb", "w")

        self.mol_one = Chem.MolFromSmiles("CCC[NH+]CCC.Cl")
        self.mol_two = Chem.MolFromSmiles("CCC[NH+]CCC")
        self.mol_three = Chem.MolFromSmiles("CCCNCCC")
        self.mol_sd_str = "\n     RDKit          3D\n\n 18 20  0  0  0  0  0  0  0  0999 V2000\n  -36.5750   18.3890   81.7690 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.7570   21.7830   81.5850 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -41.9200   21.2060   80.8120 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -41.5640   19.9910   79.9780 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.8550   18.8980   80.7650 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -32.6180   19.7650   80.0400 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -35.2150   18.4150   82.0930 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -34.3720   18.9670   81.1730 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -34.8360   19.4860   79.9810 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -36.1590   19.4760   79.6490 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -37.0540   18.9200   80.5800 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -38.5110   18.8910   80.1950 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -39.2240   19.9070   82.3780 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.2870   20.9350   82.7490 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -39.4600   19.2230   81.0980 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -38.8130   18.5540   79.0490 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -33.7840   19.9820   79.2450 O   0  0  0  0  0  0  0  0  0  0  0  0\n  -33.0050   19.1080   81.2490 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  7  2  0\n  1 11  1  0\n  2  3  1  0\n  2 14  1  0\n  3  4  1  0\n  4  5  1  0\n  5 15  1  0\n  6 17  1  0\n  6 18  1  0\n  7  8  1  0\n  8  9  2  0\n  8 18  1  0\n  9 10  1  0\n  9 17  1  0\n 10 11  2  0\n 11 12  1  0\n 12 15  1  0\n 12 16  2  0\n 13 14  1  0\n 13 15  1  0\nM  END\n"
        out_f = open("mol.sd_one", "w")
        out_f.write(self.mol_sd_str)
        self.mol_sd = out_f.name
        self.mol_sd_two_str = "\n     RDKit          3D\n\n 18 19  0  0  0  0  0  0  0  0999 V2000\n  -34.3680   23.4240   76.0460 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.7600   21.3970   82.8690 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.3490   20.2370   83.4190 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.1630   19.0970   82.7260 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.3690   19.1280   81.3680 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -35.7060   24.0440   75.8520 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -36.8250   22.2590   77.2910 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -38.0350   22.2270   77.9400 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -38.6470   23.4280   77.5950 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -38.7950   21.1690   78.6150 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -41.0140   20.3310   79.2620 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.8030   20.2790   80.7470 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.9790   21.4150   81.5160 C   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.1360   20.2240   84.7490 F   0  0  0  0  0  0  0  0  0  0  0  0\n  -36.7580   23.4430   76.6640 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -37.8780   24.1510   76.8320 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -40.1350   21.3160   78.6610 N   0  0  0  0  0  0  0  0  0  0  0  0\n  -38.2430   20.1280   78.9310 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0\n  2  3  2  0\n  2 13  1  0\n  3  4  1  0\n  3 14  1  0\n  4  5  2  0\n  5 12  1  0\n  6 15  1  0\n  7  8  2  0\n  7 15  1  0\n  8  9  1  0\n  8 10  1  0\n  9 16  2  0\n 10 17  1  0\n 10 18  2  0\n 11 12  1  0\n 11 17  1  0\n 12 13  2  0\n 15 16  1  0\nM  END\n"
        out_f = open("mol.sd_one_two", "w")
        out_f.write(self.mol_sd_two_str)
        self.mol_sd_two = out_f.name
        self.mol_smi =  "O=C(c1ccc2c(c1)OCO2)N1CCCCCC1"
        self.mol_smi_two = "Cn1cc(C(=O)NCc2ccc(F)cc2)cn1"
        # Build the dummy database file
        self.target = Target.objects.create(title="DUMMY_TARGET")
        self.cmpd = Compound.objects.create(inchi="DUM_INCH",smiles=self.mol_smi,mol_log_p=0.1,
                                            mol_wt=0.2,tpsa=0.3,heavy_atom_count=1,heavy_atom_mol_wt=2,
                                            nhoh_count=3,no_count=4,num_h_acceptors=5,num_h_donors=6,
                                            num_het_atoms=7,num_rot_bonds=8,num_val_electrons=9,ring_count=10)
        self.cmpd_two = Compound.objects.create(inchi="DUM_INCH_TWO", smiles=self.mol_smi_two, mol_log_p=0.1,
                                            mol_wt=0.2, tpsa=0.3, heavy_atom_count=1, heavy_atom_mol_wt=2,
                                            nhoh_count=3, no_count=4, num_h_acceptors=5, num_h_donors=6,
                                            num_het_atoms=7, num_rot_bonds=8, num_val_electrons=9, ring_count=10)
        self.protein = Protein.objects.create(code="DUMM",target_id=self.target,
                                              pdb_info="my_pdb.pdb")
        self.protein_two = Protein.objects.create(code="DUMM_TWO", target_id=self.target,
                                              pdb_info="my_pdb.pdb")
        self.mol = Molecule.objects.create(smiles="O=C(c1ccc2c(c1)OCO2)N1CCCCCC1",lig_id="DUM",chain_id="C",
                                           sdf_info=self.mol_sd_str,
                                rscc=0.1,occupancy=0.2,x_com=0.3,y_com=0.4,z_com=0.5,rmsd=0.6,
                                prot_id=self.protein,cmpd_id=self.cmpd)
        self.dj_mol_two = Molecule.objects.create(smiles=self.mol_smi_two, lig_id="DUM", chain_id="C",
                                           sdf_info=self.mol_sd_two_str,
                                           rscc=0.1, occupancy=0.2, x_com=0.3, y_com=0.4, z_com=0.5, rmsd=0.6,
                                           prot_id=self.protein_two, cmpd_id=self.cmpd_two)

    def test_add_target(self):
        add_target("DUMMY_TARGET")
        add_target("DUMMY_TWO")

    def test_add_mol(self):
        self.assertEqual(Molecule.objects.count(),2)
        add_mol(self.mol_sd_two,self.protein_two)
        self.assertEqual(Molecule.objects.count(),3)
        add_mol(self.mol_sd_two,self.protein_two)
        self.assertEqual(Molecule.objects.count(),3)
        add_mol(self.mol_sd,self.protein)
        self.assertEqual(Molecule.objects.count(),4)
        add_mol(self.mol_sd, self.protein)
        self.assertEqual(Molecule.objects.count(), 4)

    def test_add_comp(self):
        self.assertEqual(Compound.objects.count(),2)
        add_comp(self.mol_one)
        self.assertEqual(Compound.objects.count(),3)
        add_comp(self.mol_two)
        self.assertEqual(Compound.objects.count(),3)
        add_comp(self.mol_three)
        self.assertEqual(Compound.objects.count(),3)

    def test_add_prot(self):
        self.assertEqual(Protein.objects.count(),2)
        add_prot("DUMMY_NEW.pdb","DUMMY_NEW",self.target)
        self.assertEqual(Protein.objects.count(),3)
        self.assertTrue(Protein.objects.get(code="DUMMY_NEW").pdb_info.path.startswith("/code/media/pdbs/DUMMY_NEW"))
        add_prot("DUMMY_NEW_TWO.pdb","DUMMY_NEW",self.target)
        self.assertEqual(Protein.objects.count(),3)
        self.assertTrue(Protein.objects.get(code="DUMMY_NEW").pdb_info.path.startswith("/code/media/pdbs/DUMMY_NEW_TWO"))

    def test_analyse_mols(self):
        analyse_mols(Molecule.objects.filter(prot_id__target_id=self.target),self.target)
        self.assertEqual(MolGroup.objects.filter(group_type='MC').count(),1)
        self.assertEqual(MolGroup.objects.filter(group_type='PC').count(),16)
        analyse_mols(Molecule.objects.filter(prot_id__target_id=self.target),self.target)
        self.assertEqual(MolGroup.objects.filter(group_type='MC').count(), 2)
        self.assertEqual(MolGroup.objects.filter(group_type='PC').count(), 32)

    def test_anaylse_target(self):
        analyse_target("DUMMY_TARGET")
        self.assertEqual(MolGroup.objects.filter(group_type='MC').count(),1)
        self.assertEqual(MolGroup.objects.filter(group_type='PC').count(),16)
        analyse_target("DUMMY_TARGET")
        self.assertEqual(MolGroup.objects.filter(group_type='MC').count(),1)
        self.assertEqual(MolGroup.objects.filter(group_type='PC').count(),16)