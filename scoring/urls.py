from django.conf.urls import url

from . import views

urlpatterns = [
    url(r"^gen_conf_from_vect/$", views.gen_conf_from_vect, name="gen_conf_from_vect"),
]
