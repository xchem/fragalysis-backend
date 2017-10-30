from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^query/$', views.full_graph, name='full_graph'),
    url(r'^query/$', views.pick_mols, name='pick_mols'),
    url(r'^display/$', views.display, name='display'),

]