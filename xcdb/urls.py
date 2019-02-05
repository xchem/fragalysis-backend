from django.conf.urls import include, url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter

from xchem_db.views import TargetView, CompoundsView, ReferenceView, SoakdbFilesView, CrystalView, DataProcessingView, \
    DimpleView, LabView, RefinementView, PanddaAnalysisView, PanddaRunView, PanddaSiteView, PanddaEventView, \
    ProasisOutView

router = DefaultRouter()
# Register the basic data
router.register(r"target", TargetView)
router.register(r"compound", CompoundsView)
router.register(r'reference', ReferenceView)
router.register(r'soakdb', SoakdbFilesView)
router.register(r'crystal', CrystalView)
router.register(r'dataproc', DataProcessingView)
router.register(r'dimple', DimpleView)
router.register(r'lab', LabView)
router.register(r'refinement', RefinementView)
router.register(r'pandda_analysis', PanddaAnalysisView)
router.register(r'pandda_run', PanddaRunView)
router.register(r'pandda_site', PanddaSiteView)
router.register(r'pandda_event', PanddaEventView)
router.register(r'proasis_out', ProasisOutView)

urlpatterns = [
    url(r"^", include(router.urls)),
    url(r"^auth$", drf_views.obtain_auth_token, name="auth"),
]