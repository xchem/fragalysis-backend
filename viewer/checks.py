"""Django system checks for the viewer app.

Django runs system checks before management commands, such as the
"manage.py migrate" that the stack runs before it starts. A failed check
stops the stack early, with a clear message, instead of failing later.
"""
from botocore.exceptions import BotoCoreError, ClientError
from django.conf import settings
from django.core.checks import Error, register
from django.core.exceptions import ImproperlyConfigured
from django.core.files.storage import storages

_MEDIA_STORAGES = ("FILESYSTEM", "S3")


@register()
def check_media_storage(app_configs, **kwargs):  # pylint: disable=unused-argument
    """Checks the MEDIA_STORAGE setting and, for S3, that the bucket exists."""
    if settings.MEDIA_STORAGE not in _MEDIA_STORAGES:
        return [
            Error(
                f"MEDIA_STORAGE '{settings.MEDIA_STORAGE}' is not supported",
                hint=f"Use one of {', '.join(_MEDIA_STORAGES)}",
                id="viewer.E001",
            )
        ]
    if settings.MEDIA_STORAGE != "S3":
        return []

    try:
        storage = storages["default"]
    except ImproperlyConfigured as exc:
        return [Error(str(exc), hint="Set MEDIA_S3_BUCKET_NAME", id="viewer.E002")]

    try:
        storage.client.head_bucket(Bucket=storage.bucket_name)
    except (BotoCoreError, ClientError) as exc:
        return [
            Error(
                f"The media bucket '{storage.bucket_name}' cannot be used: {exc}",
                hint="Check that the bucket exists, and check the MEDIA_S3_*"
                " endpoint, region and credential settings",
                id="viewer.E003",
            )
        ]
    return []
