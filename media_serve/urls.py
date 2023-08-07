from django.urls import path

from . import views

urlpatterns = [
    path("pdbs/<file_path>/", views.prot_download, name="get_protein"),
    path("bound/<file_path>/", views.bound_download, name="get_bound"),
    path("metadata/<file_path>/", views.metadata_download, name="get_metadata"),
    path("targets/<file_path>/", views.archive_download, name="get_archive"),
    path("maps/<file_path>/", views.map_download, name='get_map'),
]
