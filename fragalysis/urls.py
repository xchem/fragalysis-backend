"""composeexample URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import include, url
from django.contrib import admin
import django_cas_ng.views
from graphene_django.views import GraphQLView
from django.views.generic.base import RedirectView


urlpatterns = [
    url(r"^admin/", admin.site.urls),
    url(r"^viewer/", include("viewer.urls")),
    url(r"^network/", include("network.urls")),
    url(r"^api/", include("api.urls")),
    url(r"^media_test/", include("media_serve.urls")),
    url(r"^scoring/", include("scoring.urls")),
    url(r"^graphql/", GraphQLView.as_view(graphiql=True)),
    url(r"^accounts/login/", django_cas_ng.views.login, name="cas_ng_login"),
    url(r"^accounts/logout/", django_cas_ng.views.logout, name="cas_ng_logout"),
    url(
        r"^accounts/callback$",
        django_cas_ng.views.callback,
        name="cas_ng_proxy_callback",
    ),
    url(r"^$", RedirectView.as_view(url="/accounts/login/")),
]
