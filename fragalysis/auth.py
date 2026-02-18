"""Classes to override default OIDCAuthenticationBackend (Keycloak authentication)
"""
import logging

from django.conf import settings
from mozilla_django_oidc.auth import OIDCAuthenticationBackend

logger = logging.getLogger(__name__)


class KeycloakOIDCAuthenticationBackend(OIDCAuthenticationBackend):
    def verify_claims(self, claims):
        # The default implementation insists on  an email in the token's scopes.
        # We define our own OIDC_RP_SCOPES, and need to implement our own
        # 'verify_claims()' method.
        # Implemented as part of #2054
        required_scopes: str = self.get_settings(settings.OIDC_RP_SCOPES, '')
        logger.info('OIDC_RP_SCOPES="%s" claims=%s', required_scopes, claims)

        verified: bool = True
        for scope in required_scopes.split():
            if scope not in claims:
                logger.error('Claim has no "%s" (required by OIDC_RP_SCOPES)', scope)
                verified = False

        return verified

    # Overrides Authentication Backend so that Django users are created with the keycloak preferred_username
    def create_user(self, claims):
        user = super(KeycloakOIDCAuthenticationBackend, self).create_user(claims)
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
        username = claims.get(settings.SCOPE_USERNAME_FIELD)
        assert username

        user.username = username
        user.email = claims.get('email', '')
        user.first_name = claims.get('given_name', '')
        user.last_name = claims.get('family_name', '')
        user.save()
        return user
