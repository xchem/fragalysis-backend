from django.conf.urls import include, url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter

from api import views as xchem_views
from hotspots import views as hostpot_views
from hypothesis import views as hypo_views
from scoring import views as score_views
from viewer import views as viewer_views

router = DefaultRouter()
# Register the basic data
router.register(r"molecules", viewer_views.MoleculeView)
router.register(r"compounds", viewer_views.CompoundView)
router.register(r"targets", viewer_views.TargetView)
router.register(r"proteins", viewer_views.ProteinView)
# Get the derived data
router.register(r"molimg", viewer_views.MolImageView)
router.register(r"vector", viewer_views.VectorsView)
router.register(r"graph", viewer_views.GraphView)
router.register(r"cmpdimg", viewer_views.CompoundImageView)
router.register(r"protmap", viewer_views.ProteinMapInfoView)
router.register(r"protpdb", viewer_views.ProteinPDBInfoView)
# Hotspot maps
router.register(r"hotspots", hostpot_views.HotspotView)
# Register the vectors and hypothesis
router.register(r"vectors", hypo_views.VectorView)
router.register(r"vector3ds", hypo_views.Vector3DView)
router.register(r"interactions", hypo_views.InteractionView)
router.register(r"proteinres", hypo_views.ProteinResidueView)
router.register(r"targetres", hypo_views.TargetResidueView)
router.register(r"interactionpoints", hypo_views.InteractionPointView)
# Register the  choices
router.register(r"scorechoice", score_views.ScoreChoiceView)
router.register(r"molchoice", score_views.MolChoiceView)
router.register(r"protchoice", score_views.ProtChoiceView)
router.register(r"cmpdchoice", score_views.CmpdChoiceView)
# Register the scenese
router.register(r"viewscene", score_views.ViewSceneView)
# Register the groups
router.register(r"molgroup", score_views.MolGroupView)
# Register the basic urls
router.register(r"crystal", xchem_views.CrystalView)
router.register(r"dataproc", xchem_views.DataProcessingView)
router.register(r"dimple", xchem_views.DimpleView)
router.register(r"lab", xchem_views.LabView)
router.register(r"refinement", xchem_views.RefinementView)
router.register(r"pandda_analysis", xchem_views.PanddaAnalysisView)
router.register(r"pandda_run", xchem_views.PanddaRunView)
router.register(r"pandda_site", xchem_views.PanddaSiteView)
router.register(r"pandda_event", xchem_views.PanddaEventView)
router.register(r"proasis_out", xchem_views.ProasisOutView)

from rest_framework_swagger.renderers import OpenAPIRenderer, SwaggerUIRenderer
from rest_framework.decorators import api_view, renderer_classes
from rest_framework import response, schemas


@api_view()
@renderer_classes([SwaggerUIRenderer, OpenAPIRenderer])
def schema_view(request):
    url = request.build_absolute_uri().replace("http://", "https://")
    generator = schemas.SchemaGenerator(title="Fragalysis API", url=url)
    return response.Response(generator.get_schema(request=request))


urlpatterns = [
    url(r"^", include(router.urls)),
    url(r"^auth$", drf_views.obtain_auth_token, name="auth"),
    url(r"^swagger$", schema_view),
]
