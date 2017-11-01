from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^full_graph/$', views.full_graph, name='full_graph'),
    url(r'^pick_mols/$', views.pick_mols, name='pick_mols'),
    url(r'^display/$', views.display, name='display'),
    url(r'^query_db/$', views.query_db, name='query_db'),

]