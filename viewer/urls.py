from django.urls import path, re_path

from . import views

app_name = 'viewer'

urlpatterns = [
    re_path(r"^react/*", views.react, name="react"),
    path("upload_cset/", views.UploadCSet.as_view(), name="upload_cset"),
    path(
        "validate_task/<uuid:validate_task_id>/",
        views.ValidateTaskView.as_view(),
        name="validate_task",
    ),
    path(
        "upload_task/<uuid:upload_task_id>/",
        views.UploadTaskView.as_view(),
        name="upload_task",
    ),
    path(
        "update_task/<uuid:update_task_id>/",
        views.UpdateTaskView.as_view(),
        name="update_task",
    ),
    path("img_from_smiles/", views.img_from_smiles, name="img_from_smiles"),
    path("highlight_mol_diff/", views.highlight_mol_diff, name="highlight_mol_diff"),
    path("sim_search/", views.similarity_search, name="sim_search"),
    path("open_targets/", views.get_open_targets, name="get_open_targets"),
    path("compound_set/<name>/", views.cset_download, name="compound_set"),
    path("protein_set/<name>/", views.pset_download, name="protein_set"),
    path("upload_designs/", views.DSetUploadView.as_view(), name="upload_designs"),
    path("job_access/", views.JobAccessView.as_view(), name="job_access"),
    path("task_status/<uuid:task_id>/", views.TaskStatus.as_view(), name="task_status"),
    path("service_state/", views.ServiceState.as_view(), name="service_state"),
]
