from django.conf.urls import include, url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter

from fragspect.views import FragspectView

router = DefaultRouter()

router.register(r'fragspect', FragspectView)

urlpatterns = [
    url(r"^", include(router.urls)),
    url(r"^auth$", drf_views.obtain_auth_token, name="auth"),
]