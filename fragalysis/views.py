# Classes/Methods to override default OIDC Views (Keycloak authentication)
import os
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
    """A simple endpoint that returns the content of the /code/VERSION file
    for the stack version and container environment variables for the frontend
    and backend.

    The VERSION file is adjusted during the CI process and should contain
    the TAG used to create an official build. For unofficial builds
    the version is likely to contain a CI reference.
    """
    undefined_value = "undefined"

    stack_version = Path('/code/VERSION').read_text(encoding='utf-8').strip()

    # b/e and f/e origin comes form container environment variables.
    #
    # We also need to deal with empty or unset strings
    # so the get() default does not help
    be_namespace = os.envrion.get('BE_NAMESPACE')
    if not be_namespace:
        be_namespace = undefined_value

    be_image_tag = os.envrion.get('BE_IMAGE_TAG')
    if not be_image_tag:
        be_image_tag = undefined_value

    fe_namespace = os.envrion.get('FE_NAMESPACE')
    if not fe_namespace:
        fe_namespace = undefined_value

    fe_branch = os.envrion.get('FE_BRANCH')
    if not fe_branch:
        fe_branch = undefined_value

    version_response = {'version': {'backend': f'{be_namespace}.{be_image_tag}',
                                    'frontend': f'{fe_namespace}.{fe_branch}',
                                    'stack': stack_version}}
    return JsonResponse(version_response)
