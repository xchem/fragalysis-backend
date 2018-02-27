from django.conf.urls import include,url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter
from api.views import MoleculeList,TargetList,CompoundList,ProteinList


router = DefaultRouter()
router.register(r'molecules', MoleculeList)
router.register(r'targets', TargetList)
router.register(r'compounds', CompoundList)
router.register(r'proteins', ProteinList)


# The API URLs are now determined automatically by the router.
# Additionally, we include the login URLs for the browsable API.
urlpatterns = [
    url(r'^auth$', drf_views.obtain_auth_token, name='auth'),
    url(r'^', include(router.urls))
]