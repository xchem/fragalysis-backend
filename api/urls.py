from django.conf.urls import include, url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter

# from xcdb import views as xchem_views
from hotspots import views as hostpot_views
from hypothesis import views as hypo_views
from scoring import views as score_views
from viewer import views as viewer_views
from xcdb import views as xcdb_views

router = DefaultRouter()
# Register the basic data
router.register(r"molecules", viewer_views.MoleculeView)
router.register(r"compounds", viewer_views.CompoundView)
router.register(r"targets", viewer_views.TargetView, "targets")
router.register(r"proteins", viewer_views.ProteinView)
router.register(r"projects", viewer_views.ProjectView)
router.register(r"session-projects", viewer_views.SessionProjectsView)
router.register(r"snapshots", viewer_views.SnapshotsView)
router.register(r"action-type", viewer_views.ActionTypeView)
router.register(r"session-actions", viewer_views.SessionActionsView)
router.register(r"snapshot-actions", viewer_views.SnapshotActionsView)
router.register(r"compound-identifier-types", viewer_views.CompoundIdentifierTypeView)
router.register(r"compound-identifiers", viewer_views.CompoundIdentifierView)

# Compounds sets
router.register(r"compound-sets", viewer_views.ComputedSetView)
router.register(r"compound-molecules", viewer_views.ComputedMoleculesView)
router.register(r"numerical-scores", viewer_views.NumericalScoresView)
router.register(r"text-scores", viewer_views.TextScoresView)
router.register(r"compound-scores", viewer_views.CompoundScoresView)
router.register(r"compound-mols-scores", viewer_views.ComputedMolAndScoreView)

# Get the derived data
router.register(r"molimg", viewer_views.MolImageView)
router.register(r"vector", viewer_views.VectorsView)
router.register(r"graph", viewer_views.GraphView)
router.register(r"cmpdimg", viewer_views.CompoundImageView)
router.register(r"protmap", viewer_views.ProteinMapInfoView)
router.register(r"protpdb", viewer_views.ProteinPDBInfoView)
router.register(r"protpdbbound", viewer_views.ProteinPDBBoundInfoView)

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

# Get the information
router.register(r"molannotation", score_views.MolAnnotationView)

# fragspect
router.register(r"fragspect", xcdb_views.FragspectCrystalView)

# discourse posts
router.register(r"discourse_post", viewer_views.DiscoursePostView, basename='discourse_post')

# Take a dictionary and return a csv
router.register(r"dicttocsv", viewer_views.DictToCsv, basename='dicttocsv')

# tags
router.register(r"tag_category", viewer_views.TagCategoryView, basename='tag_category')
router.register(r"molecule_tag", viewer_views.MoleculeTagView, basename='molecule_tag')
router.register(r"session_project_tag", viewer_views.SessionProjectTagView, basename='session_project_tag')
router.register(r"target_molecules", viewer_views.TargetMoleculesView, basename='target_molecules')

# Download a zip file of the requested contents
router.register(r"download_structures", viewer_views.DownloadStructures, basename='download_structures')

# Squonk Jobs
router.register(r"job_file_transfer", viewer_views.JobFileTransferView, basename='job_file_transfer')
router.register(r"job_callback", viewer_views.JobCallBackView, basename='job_callback')
router.register(r"job_config", viewer_views.JobConfigView, basename='job_config')

from rest_framework_swagger.renderers import OpenAPIRenderer, SwaggerUIRenderer
from rest_framework.decorators import api_view, renderer_classes
from rest_framework import response, schemas


@api_view()
@renderer_classes([SwaggerUIRenderer, OpenAPIRenderer])
def schema_view(request):
    url = request.build_absolute_uri()
    generator = schemas.SchemaGenerator(title="Fragalysis API")
    return response.Response(generator.get_schema(request=request))


urlpatterns = [
    url(r"^", include(router.urls)),
    url(r"^auth$", drf_views.obtain_auth_token, name="auth"),
    url(r"^swagger$", schema_view),

    url(r"job_request", viewer_views.JobRequestView.as_view(), name="job_request"),
]
