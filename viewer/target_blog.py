"""Download a target blog card's image from a URL (``image_url``).

The URL comes from a user (any member of the target's proposal), so the fetch
is deliberately narrow: http(s) only, hosts on the
``TARGET_BLOG_IMAGE_URL_ALLOWED_HOSTS`` allow-list only, no redirects (a
redirect could lead off the allow-list), a timeout, and the body is read no
further than the image size limit. The bytes are returned as an uploaded file,
so they go through exactly the same validation as an uploaded image.
"""

import logging
from io import BytesIO
from pathlib import PurePosixPath
from urllib.parse import urlsplit

import requests
from django.conf import settings
from django.core.files.uploadedfile import SimpleUploadedFile
from PIL import Image, UnidentifiedImageError

logger = logging.getLogger(__name__)

_CHUNK_BYTES = 64 * 1024


class ImageFetchError(Exception):
    """The image could not be fetched. The message is safe to show the user."""


def fetch_image(url: str) -> SimpleUploadedFile:
    """Download the image at ``url``, or raise ImageFetchError."""
    parts = urlsplit(url)
    if parts.scheme not in ("http", "https"):
        raise ImageFetchError("Only http and https image URLs can be fetched")
    allowed = settings.TARGET_BLOG_IMAGE_URL_ALLOWED_HOSTS_LIST
    if not allowed:
        raise ImageFetchError("Fetching images from a URL is disabled")
    host = (parts.hostname or "").lower()
    if host not in allowed:
        raise ImageFetchError(
            f"Images can only be fetched from: {', '.join(sorted(allowed))}"
        )

    max_bytes = settings.TARGET_BLOG_IMAGE_MAX_BYTES
    try:
        with requests.get(
            url,
            stream=True,
            allow_redirects=False,
            timeout=settings.TARGET_BLOG_IMAGE_URL_TIMEOUT_S,
        ) as response:
            if response.is_redirect:
                raise ImageFetchError(
                    "The image URL redirects elsewhere, give the final URL instead"
                )
            if response.status_code != 200:
                raise ImageFetchError(
                    f"The image URL answered with HTTP {response.status_code}"
                )
            content = bytearray()
            for chunk in response.iter_content(_CHUNK_BYTES):
                content.extend(chunk)
                if len(content) > max_bytes:
                    raise ImageFetchError(f"The image is larger than {max_bytes} bytes")
    except requests.exceptions.RequestException as exc:
        logger.warning("Could not fetch blog image %s: %s", url, exc)
        raise ImageFetchError("The image URL could not be reached") from exc

    return SimpleUploadedFile(
        _file_name(parts.path, bytes(content)),
        bytes(content),
        content_type=response.headers.get("Content-Type"),
    )


def _file_name(url_path: str, content: bytes) -> str:
    """The URL's file name, given an extension if it has none (the image
    validation rejects a name without one). Not-an-image content is left for
    that validation to reject.
    """
    name = PurePosixPath(url_path).name or "image"
    if PurePosixPath(name).suffix:
        return name
    try:
        with Image.open(BytesIO(content)) as image:
            image_format = image.format
    except (UnidentifiedImageError, OSError):
        return name
    return f"{name}.{image_format.lower()}" if image_format else name
