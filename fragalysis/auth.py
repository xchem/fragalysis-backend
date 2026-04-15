"""Classes to override default OIDCAuthenticationBackend (Keycloak authentication)
"""
import logging

from django.conf import settings
from mozilla_django_oidc.auth import OIDCAuthenticationBackend

logger = logging.getLogger(__name__)


class KeycloakOIDCAuthenticationBackend(OIDCAuthenticationBackend):
    def verify_claims(self, claims):
        # Call the super-class method
        # We don't care about the result,
        # we do our own validation.
        _ = super(KeycloakOIDCAuthenticationBackend, self).verify_claims(claims)

        # The designated username field
        # must be in the token's claims map.
        if settings.OIDC_CLAIM_USERNAME_FIELD not in claims:
            logger.info("Given claims=%s", claims)
            logger.error(
                "The '%s' field is missing from the given token's claims."
                " Without this field the login cannot be considered valid.",
                settings.OIDC_CLAIM_USERNAME_FIELD,
            )
            return False

        return True

    def create_user(self, claims):
        user = super(KeycloakOIDCAuthenticationBackend, self).create_user(claims)

        logger.debug("claims=%s", claims)

        # Get 'required' properties from the claims
        username = claims.get(settings.OIDC_CLAIM_USERNAME_FIELD)
        assert username
        user.username = username
        # Optional fields...
        user.email = claims.get('email', '')
        user.first_name = claims.get('given_name', '')
        user.last_name = claims.get('family_name', '')

        user.save()
        return user

    def filter_users_by_claims(self, claims):
        """Return all users matching the specified email.
        If nothing found matching the email, then try the username
        """
        logger.debug("claims=%s", claims)

        email = claims.get('email')
        if email:
            users = self.UserModel.objects.filter(email__iexact=email)
        else:
            username = claims.get(settings.OIDC_CLAIM_USERNAME_FIELD)
            if not username:
                return self.UserModel.objects.none()
            users = self.UserModel.objects.filter(username__iexact=username)

        return users

    def update_user(self, user, claims):
        """Update a user from a claim.
        We need the expected username field.
        """
        logger.debug("user=%s (username=%s) claims=%s", user, user.username, claims)

        username = claims.get(settings.OIDC_CLAIM_USERNAME_FIELD)
        assert username

        # Log the existing user record values we're about to change (m2ms-2116)
        logger.info(
            "[-] id=%s username=%s email=%s given_name=%s family_name=%s",
            user.id,
            user.username,
            user.email,
            user.first_name,
            user.last_name,
        )
        # And the incoming values...
        new_email = claims.get('email', '')
        new_given_name = claims.get('given_name', '')
        new_last_name = claims.get('family_name', '')
        logger.info(
            "[+] id=%s username=%s email=%s given_name=%s family_name=%s",
            user.id,
            username,
            new_email,
            new_given_name,
            new_last_name,
        )

        user.username = username
        user.email = new_email
        user.first_name = new_given_name
        user.last_name = new_last_name

        user.save()
        return user
