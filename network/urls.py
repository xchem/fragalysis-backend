from django.urls import path

from . import views

urlpatterns = [
    path(r"full_graph/", views.full_graph, name="full_graph"),
]
