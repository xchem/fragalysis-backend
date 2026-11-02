# pylint: disable=redefined-outer-name,unused-argument
"""Tests for media file serving (issue 1038, option C).

Django checks access, then ``media_response()`` gives the file to nginx with
an ``X-Accel-Redirect`` to an internal location. Django does not send the
file content itself:

- FILESYSTEM storage: ``/media-local/<name>``, which nginx aliases to
  MEDIA_ROOT.
- S3 storage: ``/media-s3/<name>``, with the presigned S3 URL in
  ``X-Media-Url``. nginx proxies that URL and does not give it to the browser.

Django answers a HEAD request itself, with the file size. A URL presigned for
GET does not work for HEAD.
"""
from urllib.parse import parse_qs, urlparse

import pytest
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
from django.http import Http404
from django.test import RequestFactory

from viewer.media_storage import media_response

_NAME = "target_loader_data/Mpro/aligned/Mpro-x0001/Mpro-x0001_apo.pdb"
_CONTENT = b"ATOM      1  N   SER A   1\n"
_DISPOSITION = "attachment;filename=Mpro-x0001_apo.pdb"


def _get(name=_NAME, method="get"):
    request = getattr(RequestFactory(), method)("/media/" + name)
    return media_response(
        request,
        name,
        content_type="application/x-pilot",
        content_disposition=_DISPOSITION,
    )


@pytest.fixture
def stored_file(media_backend):
    """Stores a file with the given name (default ``_NAME``) in the media storage."""

    def _store(name=_NAME, content=_CONTENT):
        saved = default_storage.save(name, ContentFile(content))
        assert saved == name
        return name

    return _store


# GET ------------------------------------------------------------------------


@pytest.mark.parametrize("media_backend", ["filesystem"], indirect=True)
def test_filesystem_get_redirects_to_the_local_location(media_backend):
    response = _get()

    assert response.status_code == 200
    assert response["X-Accel-Redirect"] == "/media-local/" + _NAME
    assert "X-Media-Url" not in response
    assert response["Content-Type"] == "application/x-pilot"
    assert response["Content-Disposition"] == _DISPOSITION
    assert response.content == b""


@pytest.mark.parametrize("media_backend", ["s3"], indirect=True)
def test_s3_get_redirects_to_the_s3_location_with_a_presigned_url(media_backend):
    response = _get()

    assert response.status_code == 200
    assert response["X-Accel-Redirect"] == "/media-s3/" + _NAME
    assert response["Content-Type"] == "application/x-pilot"
    assert response["Content-Disposition"] == _DISPOSITION
    assert response.content == b""

    url = urlparse(response["X-Media-Url"])
    assert url.scheme == "https"
    assert url.path.endswith("/" + _NAME)
    assert "X-Amz-Signature" in parse_qs(url.query)


@pytest.mark.parametrize("media_backend", ["s3"], indirect=True)
def test_s3_presigned_url_includes_the_storage_location(media_backend, settings):
    default = settings.STORAGES["default"]
    settings.STORAGES = {
        **settings.STORAGES,
        "default": {**default, "OPTIONS": {**default["OPTIONS"], "location": "media"}},
    }

    response = _get()

    assert urlparse(response["X-Media-Url"]).path.endswith("/media/" + _NAME)


@pytest.mark.parametrize(
    "name", ["computed_set_data/a b#1.pdb", "pdbs/x?y.pdb", "pdbs/50%.pdb"]
)
def test_unsafe_characters_in_the_name_are_quoted(media_backend, name):
    response = _get(name)

    location = "/media-local/" if media_backend == "filesystem" else "/media-s3/"
    quoted = name.replace("%", "%25").replace(" ", "%20")
    quoted = quoted.replace("#", "%23").replace("?", "%3F")
    assert response["X-Accel-Redirect"] == location + quoted


def test_get_does_not_check_that_the_file_exists(media_backend):
    """nginx reports a missing file (404). Django does not ask the storage."""
    response = _get("pdbs/missing.pdb")

    assert response.status_code == 200


# HEAD -----------------------------------------------------------------------


def test_head_is_answered_with_the_file_size(stored_file):
    stored_file()

    response = _get(method="head")

    assert response.status_code == 200
    assert response["Content-Length"] == str(len(_CONTENT))
    assert response["Content-Type"] == "application/x-pilot"
    assert response["Content-Disposition"] == _DISPOSITION
    assert "X-Accel-Redirect" not in response
    assert "X-Media-Url" not in response


def test_head_of_a_missing_file_is_not_found(media_backend):
    with pytest.raises(Http404):
        _get("pdbs/missing.pdb", method="head")


# Unsafe names ---------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["", "/etc/passwd", "../secret", "pdbs/../../secret", "pdbs/.."]
)
@pytest.mark.parametrize("method", ["get", "head"])
def test_names_outside_the_media_storage_are_not_found(media_backend, name, method):
    with pytest.raises(Http404):
        _get(name, method=method)
