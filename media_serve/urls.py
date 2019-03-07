from django.conf.urls import url

from . import views

urlpatterns = [url(r"^pdbs/(?P<file_path>.+)", views.prot_download, name="get_protein"),
               url(r"^bound/(?P<file_path>.+)", views.bound_download, name="get_bound"),
               url(r"^maps/(?P<file_path>.+)", views.map_download, name='get_map')]
