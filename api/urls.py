from django.conf.urls import include,url
from rest_framework.authtoken import views as drf_views
from rest_framework.routers import DefaultRouter
from api import views


urlpatterns = [
    url(r'^auth$', drf_views.obtain_auth_token, name='auth'),
    url(r'^molecules/(?P<pk>[0-9]+)/$', views.MoleculeList.as_view()),
    url(r'^compounds/(?P<pk>[0-9]+)/$', views.CompoundList.as_view()),
    url(r'^proteins/(?P<pk>[0-9]+)/$', views.ProteinList.as_view()),
    url(r'^targets/(?P<pk>[0-9]+)/$', views.TargetList.as_view()),
]