"""A Django storage class for S3 (or any S3-compatible store), built on boto3.

The class is configured only through the ``OPTIONS`` of a Django ``STORAGES``
entry, for example::

    STORAGES = {
        "default": {
            "BACKEND": "frag_s3_storage.S3Storage",
            "OPTIONS": {"bucket_name": "fragalysis-media", "location": "media"},
        },
    }

Credentials that are not given as options come from boto3's normal credential
chain (environment variables, a shared credentials file, an instance role).
Do not pass credentials to an ``S3Storage`` used directly in a model field:
the field's migration would then record them.

This module must not import from the backend's apps. It depends only on
Django and boto3, so that it can move to its own repository unchanged.
"""
import mimetypes
import posixpath
import tempfile
from typing import IO, Any

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured, SuspiciousFileOperation
from django.core.files.base import File
from django.core.files.storage import Storage
from django.utils import timezone
from django.utils.deconstruct import deconstructible

# The S3 error codes for an object that does not exist.
# head_object() reports "404", get_object() reports "NoSuchKey".
_NOT_FOUND_CODES = {"404", "NoSuchKey"}

# Content larger than this goes to a disk file, not memory.
_SPOOL_MAX_BYTES = 10 * 1024 * 1024


def _is_not_found(error: ClientError) -> bool:
    return error.response.get("Error", {}).get("Code") in _NOT_FOUND_CODES


@deconstructible(path="frag_s3_storage.S3Storage")
class S3Storage(Storage):
    """Stores files as objects in an S3 bucket.

    S3 has no real directories. A key such as ``a/b/c.txt`` is one object,
    and ``listdir()`` derives "directories" from the ``/`` separators.
    """

    def __init__(
        self,
        bucket_name: str | None = None,
        location: str = "",
        endpoint_url: str | None = None,
        region_name: str | None = None,
        access_key_id: str | None = None,
        secret_access_key: str | None = None,
        url_expiry: int = 3600,
    ) -> None:
        if not bucket_name:
            raise ImproperlyConfigured("S3Storage needs a 'bucket_name' option")
        self.bucket_name = bucket_name
        self.location = location.strip("/")
        self.url_expiry = url_expiry
        self._client_kwargs = {
            "endpoint_url": endpoint_url,
            "region_name": region_name,
            "aws_access_key_id": access_key_id,
            "aws_secret_access_key": secret_access_key,
            "config": Config(signature_version="s3v4"),
        }
        self._client: Any = None

    @property
    def client(self) -> Any:
        """The boto3 S3 client, created on first use."""
        if self._client is None:
            self._client = boto3.client("s3", **self._client_kwargs)
        return self._client

    def _key(self, name: str) -> str:
        """Returns the object key for a storage name.

        Names that resolve outside the storage location are rejected.
        """
        name = name.replace("\\", "/")
        clean = posixpath.normpath(name) if name else "."
        if clean == ".":
            return self.location
        if clean.startswith("/") or clean == ".." or clean.startswith("../"):
            raise SuspiciousFileOperation(
                f"The name '{name}' is outside the storage location"
            )
        return posixpath.join(self.location, clean) if self.location else clean

    def _head(self, name: str) -> dict[str, Any]:
        try:
            return self.client.head_object(Bucket=self.bucket_name, Key=self._key(name))
        except ClientError as error:
            if _is_not_found(error):
                raise FileNotFoundError(name) from error
            raise

    def _open(self, name: str, mode: str = "rb") -> File:
        """Downloads the object and returns it as a (binary) read-only File."""
        if any(c in mode for c in "wax+"):
            raise ValueError(f"S3Storage can only open files for reading, not '{mode}'")
        spool = tempfile.SpooledTemporaryFile(max_size=_SPOOL_MAX_BYTES)
        try:
            self.client.download_fileobj(self.bucket_name, self._key(name), spool)
        except ClientError as error:
            spool.close()
            if _is_not_found(error):
                raise FileNotFoundError(name) from error
            raise
        spool.seek(0)
        return File(spool, name=name)

    def _save(self, name: str, content: File) -> str:
        """Uploads the content. Text content is stored as UTF-8."""
        if hasattr(content, "seek"):
            content.seek(0)
        content_type = mimetypes.guess_type(name)[0] or "application/octet-stream"
        with tempfile.SpooledTemporaryFile(max_size=_SPOOL_MAX_BYTES) as spool:
            self._copy_as_bytes(content, spool)
            spool.seek(0)
            self.client.upload_fileobj(
                spool,
                self.bucket_name,
                self._key(name),
                ExtraArgs={"ContentType": content_type},
            )
        return name

    @staticmethod
    def _copy_as_bytes(content: File, target: IO[bytes]) -> None:
        for chunk in content.chunks():
            target.write(chunk.encode("utf-8") if isinstance(chunk, str) else chunk)

    def delete(self, name: str) -> None:
        # S3 does not report an error for a missing key,
        # which matches Django's contract for delete().
        self.client.delete_object(Bucket=self.bucket_name, Key=self._key(name))

    def exists(self, name: str) -> bool:
        try:
            self._head(name)
        except FileNotFoundError:
            return False
        return True

    def size(self, name: str) -> int:
        return self._head(name)["ContentLength"]

    def get_modified_time(self, name: str):
        modified = self._head(name)["LastModified"]
        return modified if settings.USE_TZ else timezone.make_naive(modified)

    def listdir(self, path: str) -> tuple[list[str], list[str]]:
        prefix = self._key(path)
        if prefix:
            prefix += "/"
        directories: list[str] = []
        files: list[str] = []
        paginator = self.client.get_paginator("list_objects_v2")
        for page in paginator.paginate(
            Bucket=self.bucket_name, Prefix=prefix, Delimiter="/"
        ):
            for common in page.get("CommonPrefixes", []):
                directories.append(common["Prefix"][len(prefix) :].rstrip("/"))
            for obj in page.get("Contents", []):
                files.append(obj["Key"][len(prefix) :])
        return directories, files

    def url(self, name: str) -> str:
        """Returns a presigned GET URL that expires after ``url_expiry`` seconds."""
        return self.client.generate_presigned_url(
            "get_object",
            Params={"Bucket": self.bucket_name, "Key": self._key(name)},
            ExpiresIn=self.url_expiry,
        )
