"""Access to the backend's media files, for either media storage (issue 1038).

MEDIA_STORAGE selects the storage: the local filesystem (under MEDIA_ROOT) or
an S3 bucket (see ``frag_s3_storage``). Code that serves media files uses this
module, so that it does not need to know which storage is in use.

Media files are served by nginx, not by Django. A view checks that the user
can see the file, then returns ``media_response()``. That response has an
``X-Accel-Redirect`` to an internal nginx location (see "django_nginx.conf"):

- ``/media-local/<name>``: nginx reads the file from MEDIA_ROOT.
- ``/media-s3/<name>``: nginx proxies the presigned S3 URL in the
  ``X-Media-Url`` header. The browser never sees that URL, so the file stays
  on the same origin (the page is cross-origin isolated) and the URL
  cannot be shared.
"""
import posixpath
from urllib.parse import quote

from django.core.files.storage import storages
from django.http import Http404, HttpRequest, HttpResponse

from frag_s3_storage import S3Storage

# The internal nginx locations, see "django_nginx.conf".
_LOCAL_LOCATION = "/media-local/"
_S3_LOCATION = "/media-s3/"


def _safe_name(name: str) -> str:
    """Returns the name if it is inside the media storage, otherwise raises Http404."""
    clean = posixpath.normpath(name) if name else ""
    if (
        not clean
        or clean in (".", "..")
        or clean.startswith(("/", "../"))
        or "\\" in clean
    ):
        raise Http404(f"'{name}' is not a media file")
    return clean


def media_response(
    request: HttpRequest,
    name: str,
    *,
    content_type: str,
    content_disposition: str,
) -> HttpResponse:
    """Returns the response that serves a media file.

    ``name`` is the file's storage name (i.e. a FileField's ``name``).
    Call this only after the user's access to the file is checked.

    nginx serves a GET request. Django answers a HEAD request itself
    (with the file size), because a URL presigned for GET does not work
    for HEAD. A HEAD request for a missing file raises Http404.
    """
    name = _safe_name(name)
    storage = storages["default"]

    response = HttpResponse()
    response["Content-Type"] = content_type
    response["Content-Disposition"] = content_disposition

    if request.method == "HEAD":
        try:
            response["Content-Length"] = str(storage.size(name))
        except FileNotFoundError as exc:
            raise Http404(f"'{name}' does not exist") from exc
        return response

    # Quote the name: a '?' or '#' would otherwise end the nginx path.
    quoted_name = quote(name)
    if isinstance(storage, S3Storage):
        response["X-Accel-Redirect"] = _S3_LOCATION + quoted_name
        response["X-Media-Url"] = storage.presigned_url(name)
    else:
        response["X-Accel-Redirect"] = _LOCAL_LOCATION + quoted_name
    return response
