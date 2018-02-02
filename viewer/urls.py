from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^display/$', views.display, name='display'),
    url(r'^mol_view/$', views.mol_view, name='mol_view'),
    url(r'^post_view/$', views.post_view, name='post_view'),
]