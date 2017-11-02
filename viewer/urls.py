from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^display/$', views.display, name='display'),
]