from rest_framework.test import APIRequestFactory
from rest_framework.test import force_authenticate
from django.test import TestCase
from api.utils import draw_mol,_transparentsvg,get_token
# Test all these functions



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