from django.conf.urls import url

from . import views

urlpatterns = [
    url(r"^react/*", views.react, name="react"),
    url(r"^upload_cset/", views.UploadCSet.as_view(), name="upload_cset"),
    url(r"^cset_key/", views.cset_key, name="cset_key"),
    url(r"^task/<str:task_id>/", views.TaskView.as_view(), name='task'),
    url(r"^img_from_smiles/$", views.img_from_smiles, name="img_from_smiles"),
    url(r"^highlight_mol_diff/$", views.highlight_mol_diff, name="highlight_mol_diff"),
    url(r"^sim_search/$", views.similarity_search, name="sim_search"),
    url(r"^open_targets/", views.get_open_targets, name="get_open_targets"),
    url(r'^compound_set/(?P<name>.+)/$', views.cset_download, name='compound_set'),
    url(r'^protein_set/(?P<name>.+)/$', views.pset_download, name='protein_set'),
]
