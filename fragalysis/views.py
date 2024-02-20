# Classes/Methods to override default OIDC Views (Keycloak authentication)
from django.conf import settings
from django.http import JsonResponse
from mozilla_django_oidc.views import OIDCLogoutView


def keycloak_logout(request):
    """Ths method is used to retrieve logout endpoint to also end the keycloak session as well as the Django session."""
    logout_endpoint = settings.OIDC_OP_LOGOUT_ENDPOINT
    return (
        logout_endpoint
        + "?redirect_uri="
        + request.build_absolute_uri(settings.LOGOUT_REDIRECT_URL)
    )


class LogoutView(OIDCLogoutView):
    """Extend standard logout view to include get method (called from URL)"""

    def get(self, request):
        return self.post(request)


def version(request):
    """A simple endpoint that returns the content of various environment variables
    used to define the origin of the code during build-time.

    The VERSION file is adjusted during the CI process and should contain
    the TAG used to create an official build. For unofficial builds
    the version is likely to contain a CI reference.
    """
    # Unused args
    del request

    # b/e, f/e and stack origin comes form container environment variables.
    version_response = {
        "version": {
            "backend": f"{settings.BE_NAMESPACE}:{settings.BE_IMAGE_TAG}",
            "frontend": f"{settings.FE_NAMESPACE}:{settings.FE_IMAGE_TAG}",
            "stack": f"{settings.STACK_NAMESPACE}:{settings.STACK_VERSION}",
        }
    }
    return JsonResponse(version_response)
