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
import mozilla_django_oidc.views
from django.conf.urls import include
from django.contrib import admin
from django.urls import path, re_path
from django.views.generic.base import RedirectView
from graphene_django.views import GraphQLView

import fragalysis.views

urlpatterns = [
    re_path(r"^$", RedirectView.as_view(url="/viewer/react/landing")),
    path("version/", fragalysis.views.version, name="version"),
    path("admin/", admin.site.urls),
    path("viewer/", include("viewer.urls")),
    path("network/", include("network.urls")),
    path("api/", include("api.urls")),
    path("media/", include("media_serve.urls")),
    path("scoring/", include("scoring.urls")),
    path("xcdb/", include("xcdb.urls")),
    path("graphql/", GraphQLView.as_view(graphiql=True)),
    path('oidc/', include('mozilla_django_oidc.urls')),
    path(
        "accounts/login/",
        mozilla_django_oidc.views.OIDCAuthenticationRequestView.as_view(),
        name="keylcoak_login",
    ),
    path(
        "accounts/logout/",
        fragalysis.views.LogoutView.as_view(),
        name="keycloak_logout",
    ),
    path(
        "oidc/callback/",
        mozilla_django_oidc.views.OIDCAuthenticationCallbackView.as_view(),
        name="keycloak_callback",
    ),
    path("", include("django_prometheus.urls")),
]
