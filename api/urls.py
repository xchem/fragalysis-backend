from django.conf.urls import include,url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter
from api import views
from scoring import views as score_views
from viewer import views as viewer_views
from pandda import views as pandda_views
from hypothesis import views as hypo_views

router = DefaultRouter()
# Register the basic data
router.register(r'molecules', viewer_views.MoleculeView)
router.register(r'compounds', viewer_views.CompoundView)
router.register(r'targets', viewer_views.TargetView)
router.register(r'proteins', viewer_views.ProteinView)
router.register(r'events', pandda_views.PanddaEventView)
router.register(r'sites', pandda_views.PanddaSiteView)
# Register the vectors and hypothesis
router.register(r'vectors', hypo_views.VectorView)
router.register(r'vector3ds', hypo_views.Vector3DView)
router.register(r'interactions', hypo_views.InteractionView)
router.register(r'proteinres', hypo_views.ProteinResidueView)
router.register(r'targetres', hypo_views.TargetResidueView)
router.register(r'interactionpoint', hypo_views.InteractionPointView)
# Register the  choices
router.register(r'scorechoice',score_views.ScoreChoiceView)
router.register(r'molchoice',score_views.MolChoiceView)
router.register(r'protchoice',score_views.ProtChoiceView)
router.register(r'cmpdchoice',score_views.CmpdChoiceView)
# Register the scenese
router.register(r'viewscene',score_views.ViewSceneView)
# Register the groups
router.register(r'molgroup',score_views.MolGroupView)

urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^auth$', drf_views.obtain_auth_token, name='auth'),
]