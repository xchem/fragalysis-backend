from django.test import TestCase
from viewer.models import Molecule

class MoleculeTestCase(TestCase):

    def setUp(self):
        Molecule.objects.create()