"""Tests for the upload target-access authorisation helpers in ``api.security``.

``check_upload_tas_authorisation`` decides whether a user may upload data
against a target-access string. It is HTTP-free: it returns ``None`` when the
upload may proceed, or an ``UploadTASAuthorisationFailure`` describing why not
(the *view* turns that into a login redirect or a 403 Response).

``user_has_loader_role`` is the role predicate behind the Loader bypass shared
by the upload endpoints and ``TaskStatusView``.
"""
# pylint: disable=redefined-outer-name,unused-argument
import pytest

from api.security import (
    UploadTASAuthorisationFailure,
    check_upload_tas_authorisation,
    user_has_loader_role,
)
from viewer.models import UserRole

TAS = "lb00000-1"


@pytest.fixture
def upload_settings(settings):
    """Force the authenticated-upload path with Django-only membership."""
    settings.AUTHENTICATE_UPLOAD = True
    settings.TA_AUTH_SERVICE = ""
    return settings


@pytest.fixture
def make_request(rf):
    """A POST request carrying ``user`` and optional extra headers."""

    def _make(user, **headers):
        request = rf.post(
            "/",
            **{f"HTTP_{k.upper().replace('-', '_')}": v for k, v in headers.items()},
        )
        request.user = user
        return request

    return _make


@pytest.fixture
def make_loader(make_user):
    """A user holding the Loader role."""

    def _make(username="loader"):
        user = make_user(username)
        role, _ = UserRole.objects.get_or_create(name=UserRole.LOADER_ROLE)
        role.users.add(user)
        return user

    return _make


@pytest.mark.django_db
def test_user_has_loader_role(make_user, make_loader):
    assert user_has_loader_role(make_loader()) is True
    assert user_has_loader_role(make_user("plain")) is False


@pytest.mark.django_db
def test_authenticate_upload_disabled_allows_anyone(
    upload_settings, make_request, make_user
):
    upload_settings.AUTHENTICATE_UPLOAD = False
    request = make_request(make_user("anyone"))

    assert check_upload_tas_authorisation(request, TAS) is None


@pytest.mark.django_db
def test_unauthenticated_user_needs_login(upload_settings, make_request):
    from django.contrib.auth.models import AnonymousUser

    failure = check_upload_tas_authorisation(make_request(AnonymousUser()), TAS)

    assert isinstance(failure, UploadTASAuthorisationFailure)
    assert failure.login_required is True


@pytest.mark.django_db
def test_member_is_authorised(upload_settings, make_request, make_user, make_project):
    user = make_user("member")
    make_project(TAS, members=[user])

    assert check_upload_tas_authorisation(make_request(user), TAS) is None


@pytest.mark.django_db
def test_non_member_is_forbidden(
    upload_settings, make_request, make_user, make_project
):
    make_project(TAS)
    failure = check_upload_tas_authorisation(make_request(make_user("outsider")), TAS)

    assert failure is not None
    assert failure.login_required is False
    assert failure.error_body is not None
    assert "target_access_string" in failure.error_body


@pytest.mark.django_db
def test_loader_bypasses_membership(upload_settings, make_request, make_loader):
    """A Loader may upload to a proposal they are not a member of."""
    assert check_upload_tas_authorisation(make_request(make_loader()), TAS) is None


@pytest.mark.django_db
def test_service_account_without_django_user_header_forbidden(
    upload_settings, make_request, make_user
):
    request = make_request(make_user("asap-service"))

    failure = check_upload_tas_authorisation(request, TAS)

    assert failure is not None
    assert failure.error_body is not None
    assert "error" in failure.error_body


@pytest.mark.django_db
def test_service_account_with_unknown_django_user_forbidden(
    upload_settings, make_request, make_user
):
    request = make_request(make_user("asap-service"), **{"django-user": "nobody"})

    failure = check_upload_tas_authorisation(request, TAS)

    assert failure is not None
    assert failure.error_body is not None
    assert "error" in failure.error_body


@pytest.mark.django_db
def test_service_account_delegates_to_supplied_member(
    upload_settings, make_request, make_user, make_project
):
    member = make_user("real-member")
    make_project(TAS, members=[member])
    request = make_request(make_user("asap-service"), **{"django-user": "real-member"})

    assert check_upload_tas_authorisation(request, TAS) is None
