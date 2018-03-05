from rest_framework.test import APIRequestFactory
from rest_framework.test import force_authenticate
from django.contrib.auth.models import AnonymousUser, User
from django.test import TestCase, RequestFactory
from api.utils import draw_mol,get_token
# Test all these functions


class APIUtilesTestCase(TestCase):

    def setUp(self):
        self.factory = RequestFactory()
        self.user = User.objects.create(username="DUMMY",password="DUMMY")

    def test_can_draw(self):
        output_mol = draw_mol("C1CCCCC1")
        svg_str = '<?xml version="1.0" encoding="UTF-8"?><svg baseProfile="full" height="200px" version="1.1" width="200px" xmlns="http://www.w3.org/2000/svg" xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve">\n<rect height="200" style="opacity:1.0;fill:none;stroke:none" width="200" x="0" y="0"> </rect>\n<path d="M 190.909,100 145.455,178.73" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 190.909,100 145.455,21.2704" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 145.455,178.73 54.5455,178.73" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 54.5455,178.73 9.09091,100" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 9.09091,100 54.5455,21.2704" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n<path d="M 54.5455,21.2704 145.455,21.2704" style="fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />\n</svg>'
        self.assertEqual(output_mol,svg_str)
        self.assertTrue(type(output_mol)==str)
        none_output_mol = draw_mol("C1CcccC1")
        self.assertEqual(none_output_mol,"None Mol")


    def test_can_get_token(self):
        request = self.factory.get('/viewer/react/')
        request.user = self.user
        token_one = get_token(request)
        self.assertTrue(type(token_one)==unicode)
        self.assertNotEqual(token_one,"")
        request.user = AnonymousUser()
        token_two = get_token(request)
        self.assertNotEqual(token_two,"")
        self.assertTrue(type(token_two)==unicode)



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