from django.conf.urls import include
from django.urls import path
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
router.register("compounds", viewer_views.CompoundView, "compounds")
router.register("targets", viewer_views.TargetView, "targets")
router.register("projects", viewer_views.ProjectView)
router.register("session-projects", viewer_views.SessionProjectsView)
router.register("snapshots", viewer_views.SnapshotsView)
router.register("action-type", viewer_views.ActionTypeView)
router.register("session-actions", viewer_views.SessionActionsView)
router.register("snapshot-actions", viewer_views.SnapshotActionsView)
router.register("compound-identifier-types", viewer_views.CompoundIdentifierTypeView)
router.register("compound-identifiers", viewer_views.CompoundIdentifierView)

# Compounds sets
router.register("compound-sets", viewer_views.ComputedSetView)
router.register("compound-molecules", viewer_views.ComputedMoleculesView)
router.register("numerical-scores", viewer_views.NumericalScoresView)
router.register("text-scores", viewer_views.TextScoresView)
router.register("compound-scores", viewer_views.CompoundScoresView, "compound-scores")
router.register(
    "compound-mols-scores", viewer_views.ComputedMolAndScoreView, "compound-mols-scores"
)

# Get the derived data
router.register("molimg", viewer_views.MolImageView, "molimg")
router.register("vector", viewer_views.VectorsView, "vector")
router.register("graph", viewer_views.GraphView, "graph")
router.register("cmpdimg", viewer_views.CompoundImageView, "cmpdimg")
router.register("protmap", viewer_views.ProteinMapInfoView, "protmap")
router.register("protpdb", viewer_views.ProteinPDBInfoView, "protpdb")
router.register("protpdbbound", viewer_views.ProteinPDBBoundInfoView, "protpdbbound")
router.register("molprops", viewer_views.MolecularPropertiesView, "molprops")

# Hotspot maps
router.register("hotspots", hostpot_views.HotspotView)

# Register the vectors and hypothesis
router.register("interactions", hypo_views.InteractionView)
router.register("targetres", hypo_views.TargetResidueView)
router.register("interactionpoints", hypo_views.InteractionPointView)

# Register the  choices
router.register("scorechoice", score_views.ScoreChoiceView)
router.register("siteobservationchoice", score_views.SiteObservationChoiceView)
router.register("cmpdchoice", score_views.CmpdChoiceView)

# Register the scenese
router.register("viewscene", score_views.ViewSceneView)

# Register the groups
router.register("siteobservationgroup", score_views.SiteObservationGroupView)

# Get the information
router.register("siteobservationannotation", score_views.SiteObservationAnnotationView)

# fragspect
router.register("fragspect", xcdb_views.FragspectCrystalView)

# discourse posts
router.register(
    "discourse_post", viewer_views.DiscoursePostView, basename='discourse_post'
)

# Take a dictionary and return a csv
router.register("dicttocsv", viewer_views.DictToCsv, basename='dicttocsv')

# tags
router.register("tag_category", viewer_views.TagCategoryView, basename='tag_category')
router.register(
    "siteobservation_tag",
    viewer_views.SiteObservationTagView,
    basename='siteobservation_tag',
)
router.register(
    "session_project_tag",
    viewer_views.SessionProjectTagView,
    basename='session_project_tag',
)
router.register(
    "target_molecules", viewer_views.TargetMoleculesView, basename='target_molecules'
)

# Download a zip file of the requested contents
router.register(
    "download_structures",
    viewer_views.DownloadStructures,
    basename='download_structures',
)

# Experiments and Experiment (XChemAlign) upload support
router.register(
    "upload_target_experiments",
    viewer_views.UploadTargetExperiments,
    basename='upload_target_experiments',
)
router.register(
    "download_target_experiments",
    viewer_views.DownloadTargetExperiments,
    basename='download_target_experiments',
)


router.register(
    "target_experiment_uploads",
    viewer_views.TargetExperimentUploads,
    basename='target_experiment_uploads',
)
router.register(
    "site_observations", viewer_views.SiteObservations, basename='site_observations'
)
router.register("canon_sites", viewer_views.CanonSites, basename='canon_sites')
router.register(
    "canon_site_confs", viewer_views.CanonSiteConfs, basename='canon_site_confs'
)
router.register("xtalform_sites", viewer_views.XtalformSites, basename='xtalform_sites')

# Squonk Jobs
router.register(
    "job_file_transfer", viewer_views.JobFileTransferView, basename='job_file_transfer'
)
router.register("job_callback", viewer_views.JobCallBackView, basename='job_callback')
router.register("job_config", viewer_views.JobConfigView, basename='job_config')
router.register("job_override", viewer_views.JobOverrideView, basename='job_override')


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
    path("", include(router.urls)),
    path("auth/", drf_views.obtain_auth_token, name="auth"),
    path("swagger/", schema_view),
    path("job_request/", viewer_views.JobRequestView.as_view(), name="job_request"),
]
