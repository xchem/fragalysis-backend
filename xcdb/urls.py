from django.conf.urls import include, url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter

from xcdb.views import CrystalView, DataProcessingView, \
    DimpleView, LabView, RefinementView, PanddaEventView, \
    ProasisOutView, PanddaEventStatsView, FragspectCrystalView

router = DefaultRouter()

router.register(r'crystal', CrystalView)
router.register(r'dataproc', DataProcessingView)
router.register(r'dimple', DimpleView)
router.register(r'lab', LabView)
router.register(r'refinement', RefinementView)
router.register(r'pandda_event', PanddaEventView)
router.register(r'pandda_event_stats', PanddaEventStatsView)
router.register(r'proasis_out', ProasisOutView)
router.register(r'fragspect', FragspectCrystalView)

urlpatterns = [
    url(r"^", include(router.urls)),
    url(r"^auth$", drf_views.obtain_auth_token, name="auth"),
]