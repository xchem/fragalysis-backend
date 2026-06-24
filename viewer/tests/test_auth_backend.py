"""Tests for the Keycloak OIDC authentication backend (issue #984, phase 2).

``fragalysis.auth.KeycloakOIDCAuthenticationBackend`` carries the m2ms-2116
username-conflict recovery and the 403->PermissionDenied translation. These are
security-sensitive login paths with no previous tests, so the username
resolution, new-user creation, the 403 translation and the stale-account rename
recovery are pinned here.

The backend instantiates without contacting Keycloak (its __init__ only reads
settings), so we drive its methods directly. The single external seam - the
provider's userinfo HTTP call inside ``super().get_or_create_user`` - is patched.
"""
# pylint: disable=protected-access,redefined-outer-name
from unittest import mock

import pytest
from django.db import IntegrityError
from mozilla_django_oidc.auth import OIDCAuthenticationBackend
from requests.exceptions import HTTPError
from rest_framework.exceptions import PermissionDenied

from fragalysis.auth import (
    _CONFLICT_PREFIX,
    _DUPLICATE_KEY_MESSAGE,
    KeycloakOIDCAuthenticationBackend,
    _username_from_claims,
)


@pytest.fixture
def backend():
    return KeycloakOIDCAuthenticationBackend()


def _claims(username="alice", email="alice@example.com", **extra):
    """Build a claims dict keyed on the configured username field ('fedid')."""
    claims = {"fedid": username, "email": email}
    claims.update(extra)
    return claims


def test_username_from_primary_claim(settings):
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    assert _username_from_claims({"fedid": "alice"}) == "alice"


def test_username_falls_back_to_preferred_username(settings):
    """When the configured field is absent, 'preferred_username' is used."""
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    assert _username_from_claims({"preferred_username": "bob"}) == "bob"


def test_username_none_when_no_claim(settings):
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    assert _username_from_claims({"email": "x@example.com"}) is None


def test_verify_claims_rejects_when_no_username(backend, settings):
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    assert backend.verify_claims({"email": "x@example.com"}) is False
    assert backend.verify_claims({"fedid": "alice", "email": "x@example.com"}) is True


@pytest.mark.django_db
def test_create_user_sets_username_and_profile_from_claims(backend, settings):
    """A new login is created with the claim username, email and names."""
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    claims = _claims(
        username="alice",
        email="alice@example.com",
        given_name="Alice",
        family_name="Smith",
    )

    user = backend.create_user(claims)

    user.refresh_from_db()
    assert user.username == "alice"
    assert user.email == "alice@example.com"
    assert user.first_name == "Alice"
    assert user.last_name == "Smith"


def test_get_or_create_user_translates_403_to_permission_denied(backend, monkeypatch):
    """A 403 from the OIDC provider becomes a clean DRF PermissionDenied (not 500)."""
    response = mock.Mock(status_code=403)

    def _raise_403(*_args, **_kwargs):
        raise HTTPError(response=response)

    monkeypatch.setattr(OIDCAuthenticationBackend, "get_or_create_user", _raise_403)

    with pytest.raises(PermissionDenied):
        backend.get_or_create_user("access", "id", {})


def test_get_or_create_user_reraises_other_http_errors(backend, monkeypatch):
    """A non-403 HTTPError (e.g. 401) is re-raised for the library to handle."""
    response = mock.Mock(status_code=401)

    def _raise_401(*_args, **_kwargs):
        raise HTTPError(response=response)

    monkeypatch.setattr(OIDCAuthenticationBackend, "get_or_create_user", _raise_401)

    with pytest.raises(HTTPError):
        backend.get_or_create_user("access", "id", {})


@pytest.mark.django_db
def test_update_user_happy_path(backend, make_user, settings):
    """With no username clash the user is updated in place."""
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    user = make_user("alice", email="old@example.com")

    result = backend.update_user(
        user, _claims(username="alice", email="new@example.com", given_name="Alice")
    )

    result.refresh_from_db()
    assert result.username == "alice"
    assert result.email == "new@example.com"
    assert result.first_name == "Alice"


@pytest.mark.django_db
def test_update_user_frees_username_from_stale_account(backend, make_user, settings):
    """An email-matched login claiming a taken username renames the stale row.

    The user keeps the username they expect; the colliding (stale) account is
    moved aside with the _CONFLICT_PREFIX so the login still succeeds.
    """
    settings.OIDC_CLAIM_USERNAME_FIELD = "fedid"
    stale = make_user("bob", email="stale@example.com")
    # The real (email-matched) account that should own 'bob'.
    incoming = make_user("bob-temp", email="bob@example.com")

    result = backend.update_user(
        incoming, _claims(username="bob", email="bob@example.com")
    )

    result.refresh_from_db()
    stale.refresh_from_db()
    assert result.pk == incoming.pk
    assert result.username == "bob"
    assert stale.username == f"{_CONFLICT_PREFIX}bob"


@pytest.mark.django_db
def test_save_conflict_unrecoverable_integrity_error_becomes_403(
    backend, make_user, monkeypatch
):
    """A non-collision IntegrityError fails the login as 403, not 500."""
    user = make_user("carol")

    def _raise_other(*_args, **_kwargs):
        raise IntegrityError("null value in column violates not-null constraint")

    monkeypatch.setattr(user, "save", _raise_other)

    with pytest.raises(PermissionDenied):
        backend._save_resolving_username_conflict(user, "carol")


@pytest.mark.django_db
def test_save_conflict_without_conflicting_row_becomes_403(
    backend, make_user, monkeypatch
):
    """The collision text is present but no other row holds the name: 403, not 500."""
    user = make_user("dave")

    def _raise_duplicate(*_args, **_kwargs):
        raise IntegrityError(_DUPLICATE_KEY_MESSAGE)

    monkeypatch.setattr(user, "save", _raise_duplicate)

    # No other User owns 'unique-name', so recovery is impossible.
    with pytest.raises(PermissionDenied):
        backend._save_resolving_username_conflict(user, "unique-name")
