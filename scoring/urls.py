from django.urls import path

from . import views

urlpatterns = [
    path("gen_conf_from_vect/", views.gen_conf_from_vect, name="gen_conf_from_vect"),
    path("get_current_user_id/", views.get_current_user_id, name="get_current_user_id"),
]
