from django.conf.urls import include,url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter
from api import views
from scoring import views as score_views

router = DefaultRouter()
# Register the basic data
router.register(r'molecules', views.MoleculeView)
router.register(r'compounds', views.CompoundView)
router.register(r'targets', views.TargetView)
router.register(r'proteins', views.ProteinView)
router.register(r'events', views.PanddaEventView)
router.register(r'sites', views.PanddaSiteView)
# Register the vectors and interactions
router.register(r'vectors', views.VectorView)
router.register(r'vector3ds', views.Vector3DView)
router.register(r'interactions', views.InteractionView)
router.register(r'proteinres', views.ProteinResidueView)
router.register(r'targetres', views.TargetResidueView)
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