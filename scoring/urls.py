from django.conf.urls import url

from . import views

urlpatterns = [
    url(r"^crystal/$", views.CrystalView, name="crystal"),
    url(r"^data_processing/$", views.DataProcessingView, name="data_processing"),
    url(r"^dimple/$", views.DimpleView, name="dimple"),
]


# # Register the basic urls
# router.register(r"crystal", xchem_views.CrystalView)
# router.register(r"data_processing", xchem_views.DataProcessingView)
# router.register(r"dimple", xchem_views.DimpleView)
# router.register(r"lab", xchem_views.LabView)
# router.register(r"refinement", xchem_views.RefinementView)
# router.register(r"pandda_analysis", xchem_views.PanddaAnalysisView)
# router.register(r"pandda_run", xchem_views.PanddaRunView)
# router.register(r"pandda_site", xchem_views.PanddaSiteView)
# router.register(r"pandda_event", xchem_views.PanddaEventView)
# router.register(r"proasis_out", xchem_views.ProasisOutView)
