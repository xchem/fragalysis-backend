from django.conf.urls import url

from . import views

urlpatterns = [url(r"^pdbs/(?P<file_path>.+)", views.prot_download, name="get_protein")]
