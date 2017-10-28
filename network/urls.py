from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^query/$', views.query, name='query'),
    url(r'^display$', views.display, name='display'),

]