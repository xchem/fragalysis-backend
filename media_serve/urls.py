from django.urls import path, re_path

from . import views

urlpatterns = [
    # re_path(r"^pdbs/(?P<file_path>.+)", views.prot_download, name="get_protein"),
    path("bound/<file_path>/", views.bound_download, name="get_bound"),
    path("metadata/<file_path>/", views.metadata_download, name="get_metadata"),
    path("targets/<file_path>/", views.archive_download, name="get_archive"),
    path("maps/<file_path>/", views.map_download, name='get_map'),
    re_path(
        r"^target_loader_data/(?P<file_path>.+)", views.tld_download, name="get_tld"
    ),
    re_path(
        r"^computed_set_data/(?P<file_path>.+)", views.cspdb_download, name="get_cspdb"
    ),
    re_path(r"^pdbs/(?P<file_path>.+)", views.file_download, name="get_file"),
]
