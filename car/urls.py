from django.conf.urls import url
from rest_framework import routers

# Import standard views
from .api import (
    ProjectViewSet,
    MculeQuoteViewSet,
    BatchViewSet,
    TargetViewSet,
    MethodViewSet,
    ReactionViewSet,
    ProductViewSet,
)

# Import action views
from .api import (
    AnalyseActionViewSet,
    AddActionViewSet,
    ExtractActionViewSet,
    FilterActionViewSet,
    QuenchActionViewSet,
    SetTemperatureActionViewSet,
    StirActionViewSet,
   
)

# Import OT Session views
from .api import (
    OTSessionViewSet,
    DeckViewSet,
    PipetteViewSet,
    TipRackViewSet,
    PlateViewSet,
    WellViewSet,
    CompoundOrderViewSet,
    OTScriptViewSet,
)

from .views import UploadProject, ValidateTaskView, UploadTaskView

# Register standard routes
router = routers.DefaultRouter()
router.register("api/projects", ProjectViewSet, "projects")
router.register("api/mculequotes", MculeQuoteViewSet, "mculequotes")
router.register("api/batches", BatchViewSet, "batches")
router.register("api/targets", TargetViewSet, "targets")
router.register("api/methods", MethodViewSet, "methods")
router.register("api/reactions", ReactionViewSet, "reactions")
router.register("api/products", ProductViewSet, "products")

# Register action routes
router.register("api/analyseactions", AnalyseActionViewSet, "analyseactions")

router.register("api/addactions", AddActionViewSet, "addactions")
router.register("api/extractactions", ExtractActionViewSet, "extractactions")
router.register("api/filteractions", FilterActionViewSet, "filteractions")
router.register("api/quenchactions", QuenchActionViewSet, "quenchactions")
router.register(
    "api/set-temperatureactions",
    SetTemperatureActionViewSet,
    "set-temperatureactions",
)
router.register("api/stiractions", StirActionViewSet, "stiractions")

# Register Ot Session routes
router.register("api/otsessions", OTSessionViewSet, "otsessions")
router.register("api/decks", DeckViewSet, "decks")
router.register("api/pipettes", PipetteViewSet, "pipettes")
router.register("api/tipracks", TipRackViewSet, "tipracks")
router.register("api/plates", PlateViewSet, "plates")
router.register("api/wells", WellViewSet, "wells")
router.register("api/compoundorders", CompoundOrderViewSet, "compoundorders")
router.register("api/otscripts", OTScriptViewSet, "otscripts")

urlpatterns = [
    url("upload/", UploadProject.as_view(), name="uploadproject"),
    url(
        r"^validate_cartask/(?P<validate_task_id>.+)/$",
        ValidateTaskView.as_view(),
        name="validate_task",
    ),
    url(
        r"^upload_cartask/(?P<upload_task_id>.+)/$",
        UploadTaskView.as_view(),
        name="upload_task",
    ),
]

urlpatterns += router.urls
