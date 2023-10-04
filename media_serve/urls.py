from django.urls import path
from django.urls import re_path

from . import views

urlpatterns = [
    # path("pdbs/<file_path>/", views.prot_download, name="get_protein"),
    # path("pdbs/", views.prot_download, name="get_protein"),
    re_path(r"^pdbs/(?P<file_path>.+)", views.prot_download, name="get_protein"),
    path("bound/<file_path>/", views.bound_download, name="get_bound"),
    path("metadata/<file_path>/", views.metadata_download, name="get_metadata"),
    path("targets/<file_path>/", views.archive_download, name="get_archive"),
    path("maps/<file_path>/", views.map_download, name='get_map'),
    # path("target_loader_data/<file_path>", views.file_download, name='get_file'),
    # re_path(r"^target_loader_data/48225dbf-204a-48e1-8ae7-f1632f4dba89/Mpro-v2/Mpro/upload_2/aligned_files/Mpro_Nterm-x0029/(?P<file_path>.+)", views.file_download, name='get_file'),
    re_path(r"^target_loader_data/(?P<file_path>.+)", views.file_download, name="get_file"),
]
