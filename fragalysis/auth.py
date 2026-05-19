"""Classes to override default OIDCAuthenticationBackend (Keycloak authentication)
"""
import logging

from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.db import IntegrityError, transaction
from mozilla_django_oidc.auth import OIDCAuthenticationBackend

logger = logging.getLogger(__name__)

# Prefix applied to a pre-existing account whose username collides with the
# one an incoming (email-matched) login is claiming. Renaming the stale
# account frees the username so the user does not have to log in again.
# See m2ms-2116.
_CONFLICT_PREFIX = "m2ms-2116-"

# Only the username-collision IntegrityError carries this (Postgres) text.
# Any other IntegrityError is not something we know how to recover from.
_DUPLICATE_KEY_MESSAGE = "duplicate key value violates unique constraint"


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

    def _save_resolving_username_conflict(self, user, username):
        """Save 'user', recovering from the auth_user_username_key collision.

        'user' has been matched by email but is being assigned a username
        that another (stale) row already owns, so user.save() raises
        IntegrityError. We rename the stale row (see _CONFLICT_PREFIX) to
        free the username and retry, so the user does not have to log in
        again. The first save runs in a savepoint so the connection stays
        usable after the rollback.
        """
        try:
            with transaction.atomic():
                user.save()
            return user
        except IntegrityError as exc:
            if _DUPLICATE_KEY_MESSAGE not in str(exc):
                # A different IntegrityError (not the username collision):
                # fail the login cleanly with a 403 rather than a 500.
                logger.error(
                    "Unrecoverable IntegrityError saving id=%s username=%s: %s",
                    user.id,
                    username,
                    exc,
                )
                raise PermissionDenied(
                    "Could not complete login (database integrity error)."
                ) from exc

            conflicting = list(
                self.UserModel.objects.filter(username__iexact=username).exclude(
                    pk=user.pk
                )
            )
            if not conflicting:
                # The username collision text was present but no other row
                # holds it: we cannot recover, so 403 rather than 500.
                logger.error(
                    "Username collision for id=%s username=%s but no"
                    " conflicting row found: %s",
                    user.id,
                    username,
                    exc,
                )
                raise PermissionDenied(
                    "Could not complete login (username conflict)."
                ) from exc

            for other in conflicting:
                new_username = f"{_CONFLICT_PREFIX}{other.username}"
                # A repeat collision (prefixed name already taken) would
                # just violate the constraint again on retry, so fall
                # back to a pk-qualified name that is guaranteed unique.
                if (
                    self.UserModel.objects.filter(username__iexact=new_username)
                    .exclude(pk=other.pk)
                    .exists()
                ):
                    new_username = f"{_CONFLICT_PREFIX}{other.pk}-{other.username}"
                new_username = new_username[:150]
                logger.warning(
                    "Freeing username '%s': renaming conflicting"
                    " id=%s username=%s -> %s",
                    username,
                    other.pk,
                    other.username,
                    new_username,
                )
                other.username = new_username
                with transaction.atomic():
                    other.save()

            with transaction.atomic():
                user.save()
            return user

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

        return self._save_resolving_username_conflict(user, username)
