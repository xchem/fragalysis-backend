# Classes/Methods to override default OIDC Views (Keycloak authentication)
from pathlib import Path

from mozilla_django_oidc.views import OIDCLogoutView
from django.http import JsonResponse
from django.conf import settings


def keycloak_logout(request):
    """ Ths method is used to retrieve logout endpoint to also end the keycloak session as well as the Django session.
    """
    logout_endpoint = settings.OIDC_OP_LOGOUT_ENDPOINT
    return logout_endpoint + "?redirect_uri=" + request.build_absolute_uri(settings.LOGOUT_REDIRECT_URL)


class LogoutView(OIDCLogoutView):
    """ Extend standard logout view to include get method (called from URL)
    """

    def get(self, request):
        return self.post(request)


def version(request):
    """A simple endpoint that returns the content of the /code/VERSION file.
    The VERSION file is adjusted during the CI process and should contain
    the TAG used to create an official build. For unofficial builds
    the version is likely to contain a CI reference.
    """
    version = Path('/code/VERSION').read_text(encoding='utf-8').strip()
    return JsonResponse({'version': version})
