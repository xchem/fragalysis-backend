from django.conf.urls import url

from . import views

urlpatterns = [
    url(r"^crystal/$", views.CrystalView, name="crystal"),
    url(r"^data_processing/$", views.DataProcessingView, name="data_processing"),
    url(r"^dimple/$", views.DimpleView, name="dimple"),
]