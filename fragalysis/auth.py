"""Classes to override default OIDCAuthenticationBackend (Keycloak authentication)
"""
import logging

from django.conf import settings
from mozilla_django_oidc.auth import OIDCAuthenticationBackend

logger = logging.getLogger(__name__)


class KeycloakOIDCAuthenticationBackend(OIDCAuthenticationBackend):
    # Overrides Authentication Backend so that Django users are created with the keycloak preferred_username
    def create_user(self, claims):
        user = super(KeycloakOIDCAuthenticationBackend, self).create_user(claims)

        logger.info("claims=%s", claims)

        # Get 'expected' properties from the claims, some are optional.
        username = claims.get(settings.SCOPE_USERNAME_FIELD)
        assert username
        user.username = username
        # Optional...
        user.email = claims.get('email', '')
        user.first_name = claims.get('given_name', '')
        user.last_name = claims.get('family_name', '')
        user.save()
        return user

    def filter_users_by_claims(self, claims):
        """Return all users matching the specified email.
        If nothing found matching the email, then try the username
        """
        logger.info("claims=%s", claims)

        email = claims.get('email')
        username = claims.get(settings.SCOPE_USERNAME_FIELD)

        if email:
            users = self.UserModel.objects.filter(email__iexact=email)
        else:
            if not username:
                return self.UserModel.objects.none()
            users = self.UserModel.objects.filter(username__iexact=username)
        return users

    def update_user(self, user, claims):
        """Update a user from a claim.
        We need the expected username field.
        """
        logger.info("user=%s claims=%s", user, claims)

        username = claims.get(settings.SCOPE_USERNAME_FIELD)
        assert username

        user.username = username
        user.email = claims.get('email', '')
        user.first_name = claims.get('given_name', '')
        user.last_name = claims.get('family_name', '')
        user.save()
        return user
