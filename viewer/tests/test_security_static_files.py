"""Tests for the media-download authorisation barrier (issue #984, phase 2).

``api.security.ISPyBSafeStaticFiles`` (used by the ``/metadata/``, ``/bound/``,
``/maps/`` and ``/targets/`` downloads) is the access-control gate for static
file serving: ``get_response`` only resolves a file the requesting user is
allowed to see, otherwise it raises ``Http404``. A bug here leaks files, so the
member / non-member / public / anonymous paths are pinned here.

``ISPyBSafeStaticFiles2`` (used by ``/pdbs/``, ``/target_loader_data/`` and
``/computed_set_data/``) deliberately does *not* re-query the database - it only
builds the Nginx ``X-Accel-Redirect`` path - so its test asserts that path
construction rather than access control.

Access control is driven through the Django membership path
(``TA_AUTH_SERVICE`` unset, ``Project.user_id`` membership) - the same seam the
existing ``test_security.py`` uses.
"""
# The static-file mixins are configured by assigning attributes at call sites
# (see media_serve.views), so mypy cannot see them as declared members.
# mypy: disable-error-code="attr-defined"
# pylint: disable=redefined-outer-name
import pytest
from django.contrib.auth.models import AnonymousUser
from django.http import Http404
from django.test import RequestFactory

from api.security import ISPyBSafeStaticFiles, ISPyBSafeStaticFiles2
from viewer.models import Target

# The stored FileField value and the basename Nginx should be handed.
_METADATA_PATH = "metadata/Mpro_meta.csv"
_METADATA_NAME = "Mpro_meta.csv"


@pytest.fixture
def target_with_metadata(make_project, make_target):
    """Factory: a Target carrying a metadata file, in a (named) project.

    Returns ``(target, project)`` so the test can grant/withhold membership.
    """

    def _make(title="proposal-meta", members=(), open_to_public=False):
        project = make_project(title, members=members, open_to_public=open_to_public)
        target = make_target(project)
        target.metadata = _METADATA_PATH
        target.save()
        return target, project

    return _make


def _metadata_handler(user):
    """An ISPyBSafeStaticFiles wired up exactly as ``metadata_download`` wires it."""
    handler = ISPyBSafeStaticFiles()
    handler.model = Target
    request = RequestFactory().get("/metadata/" + _METADATA_NAME)
    request.user = user
    handler.request = request
    handler.permission_string = "project"
    handler.field_name = "metadata"
    handler.content_type = "application/x-pilot"
    handler.prefix = "/metadata/"
    handler.input_string = _METADATA_NAME
    return handler


@pytest.mark.django_db
def test_member_is_served_the_file(target_with_metadata, make_user, settings):
    """A project member gets an X-Accel-Redirect to the resolved file."""
    settings.TA_AUTH_SERVICE = ""
    user = make_user("member")
    target_with_metadata(members=[user])

    response = _metadata_handler(user).get_response()

    assert response.status_code == 200
    assert response["X-Accel-Redirect"] == "/metadata/" + _METADATA_NAME
    assert response["Content-Type"] == "application/x-pilot"
    assert _METADATA_NAME in response["Content-Disposition"]


@pytest.mark.django_db
def test_non_member_is_denied(target_with_metadata, make_user, settings):
    """A user with no membership of the target's project gets a 404, not the file."""
    settings.TA_AUTH_SERVICE = ""
    target_with_metadata(members=[])  # owned by nobody
    outsider = make_user("outsider")

    with pytest.raises(Http404):
        _metadata_handler(outsider).get_response()


@pytest.mark.django_db
def test_public_project_is_served_to_non_member(
    target_with_metadata, make_user, settings
):
    """An open_to_public project's file is served even without membership."""
    settings.TA_AUTH_SERVICE = ""
    settings.PUBLIC_TAS_LIST = []
    target_with_metadata(open_to_public=True, members=[])
    anyone = make_user("anyone")

    response = _metadata_handler(anyone).get_response()

    assert response.status_code == 200
    assert response["X-Accel-Redirect"] == "/metadata/" + _METADATA_NAME


@pytest.mark.django_db
def test_anonymous_user_denied_private_file(target_with_metadata, settings):
    """The ``user.pk is None`` path: an anonymous user gets no private file."""
    settings.TA_AUTH_SERVICE = ""
    target_with_metadata(open_to_public=False, members=[])

    with pytest.raises(Http404):
        _metadata_handler(AnonymousUser()).get_response()


@pytest.mark.django_db
def test_anonymous_user_served_public_file(target_with_metadata, settings):
    """An anonymous user may still read a public project's file."""
    settings.TA_AUTH_SERVICE = ""
    settings.PUBLIC_TAS_LIST = []
    target_with_metadata(open_to_public=True, members=[])

    response = _metadata_handler(AnonymousUser()).get_response()

    assert response.status_code == 200
    assert response["X-Accel-Redirect"] == "/metadata/" + _METADATA_NAME


def test_static_files2_builds_nginx_redirect_path():
    """ISPyBSafeStaticFiles2 maps prefix + input to an absolute X-Accel-Redirect.

    It does no DB lookup, so no access check happens here - it is the path
    builder for the loader-data downloads.
    """
    handler = ISPyBSafeStaticFiles2()
    handler.content_type = "application/x-pilot"
    handler.prefix = "/target_loader_data/"
    handler.input_string = "Mpro/aligned/Mpro-x0001/Mpro-x0001_apo.pdb"

    response = handler.get_response()

    assert response.status_code == 200
    assert (
        response["X-Accel-Redirect"]
        == "/target_loader_data/Mpro/aligned/Mpro-x0001/Mpro-x0001_apo.pdb"
    )
    assert response["Content-Type"] == "application/x-pilot"
    assert handler.input_string in response["Content-Disposition"]
