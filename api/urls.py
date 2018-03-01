from django.conf.urls import include,url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter
from api import views

router = DefaultRouter()
router.register(r'molecules', views.MoleculeView)
router.register(r'mdl', views.MDLView)
router.register(r'compounds', views.CompoundView)
router.register(r'targets', views.TargetView)
router.register(r'proteins', views.ProteinView)


urlpatterns = [
    url(r'^', include(router.urls)),
    url(r'^auth$', drf_views.obtain_auth_token, name='auth'),
]