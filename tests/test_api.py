from rest_framework.test import APIRequestFactory
from rest_framework.test import force_authenticate
from django.test import TestCase
from django.http import HttpRequest
from django.contrib.auth.models import User
from api.utils import draw_mol,get_token
# Test all these functions


class APIUtilesTestCase(TestCase):

    def setUp(self):
        user_one = User.objects.create(username="DUMMY",password="DUMMY")
        self.request_one = HttpRequest()
        self.request_one.user = user_one.get_username()

    def test_can_draw(self):
        output_mol = draw_mol("C1CCCCC1")
        self.assertEqual(output_mol,"")
        self.assertTrue(type(output_mol),str)

    def test_can_get_token(self):
        get_token(self.request_one)


# TEST ALL THESE URLS
"""
router.register(r'molecules', views.MoleculeView)
router.register(r'mdl', views.MDLView)
router.register(r'compounds', views.CompoundView)
router.register(r'targets', views.TargetView)
router.register(r'proteins', views.ProteinView)

# Register the  choices
router.register(r'scorechoice',score_views.ScoreChoiceView)
router.register(r'molchoice',score_views.MolChoiceView)
router.register(r'protchoice',score_views.ProtChoiceView)
router.register(r'cmpdchoice',score_views.CmpdChoiceView)
# Register the scenese
router.register(r'viewscene',score_views.ViewSceneView)
# Register the groups
router.register(r'molgroup',score_views.MolGroupView)

"""