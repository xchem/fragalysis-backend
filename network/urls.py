from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^full_graph/$', views.full_graph, name='full_graph'),
    url(r'^display/$', views.display, name='display'),
]