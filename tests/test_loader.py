from django.test import TestCase
from loader.functions import desalt_compound,NeutraliseCharges,get_path_or_none,sanitize_mol
from loader.loaders import add_target,add_mol,add_comp,add_prot,analyse_mols,analyse_target

# TEST ALL THE ABOVE