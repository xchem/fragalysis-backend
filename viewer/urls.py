from django.conf.urls import include, url
from rest_framework.routers import DefaultRouter

from . import views



urlpatterns = [
    url(r"^react/*", views.react, name="react"),
    url(r"^upload_cset/", views.UploadCSet.as_view(), name="upload_cset"),
    url(r"^update_cset/", views.UpdateCSet.as_view(), name="update_cset"),
    url(r"^cset_key/", views.cset_key, name="cset_key"),
    url(r"^validate_task/(?P<validate_task_id>.+)/$", views.ValidateTaskView.as_view(), name='validate_task'),
    url(r"^upload_task/(?P<upload_task_id>.+)/$", views.UploadTaskView.as_view(), name='upload_task'),
    url(r"^update_task/(?P<update_task_id>.+)/$", views.UpdateTaskView.as_view(), name='update_task'),
    url(r"^img_from_smiles/$", views.img_from_smiles, name="img_from_smiles"),
    url(r"^highlight_mol_diff/$", views.highlight_mol_diff, name="highlight_mol_diff"),
    url(r"^sim_search/$", views.similarity_search, name="sim_search"),
    url(r"^open_targets/", views.get_open_targets, name="get_open_targets"),
    url(r'^compound_set/(?P<name>.+)/$', views.cset_download, name='compound_set'),
    url(r'^protein_set/(?P<name>.+)/$', views.pset_download, name='protein_set'),
    url(r'upload_designs/', views.DSetUploadView.as_view(), name='upload_designs')
]

