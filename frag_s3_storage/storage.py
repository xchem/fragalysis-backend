"""A Django storage class for S3 (or any S3-compatible store), built on boto3.

The class is configured only through the ``OPTIONS`` of a Django ``STORAGES``
entry, for example::

    STORAGES = {
        "default": {
            "BACKEND": "frag_s3_storage.S3Storage",
            "OPTIONS": {
                "bucket_name": "fragalysis-media",
                "location": "media",
                "url_base": "/media/",
            },
        },
    }

Credentials that are not given as options come from boto3's normal credential
chain (environment variables, a shared credentials file, an instance role).
Do not pass credentials to an ``S3Storage`` used directly in a model field:
the field's migration would then record them.

This module must not import from the backend's apps. It depends only on
Django and boto3, so that it can move to its own repository unchanged.
"""
import logging
import mimetypes
import os
import posixpath
import tempfile
import time
from collections.abc import Iterable, Iterator
from datetime import datetime
from typing import IO, Any
from urllib.parse import urljoin

import boto3
from boto3.exceptions import S3UploadFailedError
from botocore.config import Config
from botocore.exceptions import ClientError
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured, SuspiciousFileOperation
from django.core.files.base import File
from django.core.files.storage import Storage
from django.utils import timezone
from django.utils.deconstruct import deconstructible
from django.utils.encoding import filepath_to_uri

logger = logging.getLogger(__name__)

# The S3 error codes for an object that does not exist.
# head_object() reports "404", get_object() reports "NoSuchKey".
_NOT_FOUND_CODES = {"404", "NoSuchKey"}

# Content larger than this goes to a disk file, not memory.
_SPOOL_MAX_BYTES = 10 * 1024 * 1024

# The largest number of keys that one DeleteObjects request accepts.
_DELETE_BATCH_SIZE = 1000


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
        url_base: str | None = None,
        max_pool_connections: int = 16,
        upload_attempts: int = 3,
        retry_delay_s: float = 1.0,
    ) -> None:
        if not bucket_name:
            raise ImproperlyConfigured("S3Storage needs a 'bucket_name' option")
        if upload_attempts < 1:
            raise ImproperlyConfigured("S3Storage needs 'upload_attempts' of 1 or more")
        # boto3 does not retry a rejected request (i.e. a 403), and an S3 store
        # can reject one request of many for no clear reason. One failed file
        # then fails a whole target upload, so upload_file() tries again.
        self.upload_attempts = upload_attempts
        self.retry_delay_s = retry_delay_s
        self.bucket_name = bucket_name
        self.location = location.strip("/")
        self.url_expiry = url_expiry
        self.url_base = url_base
        self._client_kwargs = {
            "endpoint_url": endpoint_url,
            "region_name": region_name,
            "aws_access_key_id": access_key_id,
            "aws_secret_access_key": secret_access_key,
            # A caller that transfers files with N threads needs N connections.
            "config": Config(
                signature_version="s3v4", max_pool_connections=max_pool_connections
            ),
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

    def _key_prefix(self, prefix: str) -> str:
        """Returns the key prefix for a name prefix, keeping a trailing ``/``."""
        key = self._key(prefix)
        if key and (prefix.endswith("/") or not prefix):
            key += "/"
        return key

    def _name(self, key: str) -> str:
        """Returns the storage name for an object key (the inverse of ``_key()``)."""
        return key[len(self.location) + 1 :] if self.location else key

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

    def upload_file(self, local_path: str | os.PathLike, name: str) -> str:
        """Uploads a local file to the name, and replaces any existing file.

        boto3 uses a multipart upload for a large file. A failed upload is
        tried again, up to ``upload_attempts`` times in all, with a longer
        delay each time. The last failure is raised.
        """
        content_type = mimetypes.guess_type(name)[0] or "application/octet-stream"
        extra_args = {"ContentType": content_type}
        key = self._key(name)
        for attempt in range(1, self.upload_attempts):
            try:
                self.client.upload_file(
                    str(local_path), self.bucket_name, key, ExtraArgs=extra_args
                )
                return name
            except S3UploadFailedError as error:
                logger.warning(
                    "Upload of '%s' failed (attempt %d of %d), trying again: %s",
                    name,
                    attempt,
                    self.upload_attempts,
                    error,
                )
                time.sleep(self.retry_delay_s * attempt)
        # The last attempt raises its failure.
        self.client.upload_file(
            str(local_path), self.bucket_name, key, ExtraArgs=extra_args
        )
        return name

    def download_file(self, name: str, local_path: str | os.PathLike) -> None:
        """Downloads the file to a local file (boto3 uses ranged, parallel GETs)."""
        try:
            self.client.download_file(
                self.bucket_name, self._key(name), str(local_path)
            )
        except ClientError as error:
            if _is_not_found(error):
                raise FileNotFoundError(name) from error
            raise

    def copy(self, source_name: str, target_name: str) -> None:
        """Copies a file inside the bucket (a managed copy, for any size)."""
        try:
            self.client.copy(
                {"Bucket": self.bucket_name, "Key": self._key(source_name)},
                self.bucket_name,
                self._key(target_name),
            )
        except ClientError as error:
            if _is_not_found(error):
                raise FileNotFoundError(source_name) from error
            raise

    def delete_many(self, names: Iterable[str]) -> list[str]:
        """Deletes the files, with one request for each 1000 names.

        Returns the names that S3 did not delete. A missing file is not a failure.
        """
        keys = [self._key(name) for name in names]
        failed: list[str] = []
        for start in range(0, len(keys), _DELETE_BATCH_SIZE):
            batch = keys[start : start + _DELETE_BATCH_SIZE]
            response = self.client.delete_objects(
                Bucket=self.bucket_name,
                Delete={"Objects": [{"Key": key} for key in batch], "Quiet": True},
            )
            failed.extend(self._name(e["Key"]) for e in response.get("Errors", []))
        return failed

    def iter_objects(self, prefix: str = "") -> Iterator[tuple[str, int, datetime]]:
        """Yields ``(name, size, modified)`` for every file whose name starts with
        the prefix, at any depth. Use a trailing ``/`` to list one "directory"
        (``"Mpro/"`` does not match ``"Mpro_2/a.txt"``).
        """
        paginator = self.client.get_paginator("list_objects_v2")
        for page in paginator.paginate(
            Bucket=self.bucket_name, Prefix=self._key_prefix(prefix)
        ):
            for obj in page.get("Contents", []):
                yield self._name(obj["Key"]), obj["Size"], obj["LastModified"]

    def url(self, name: str) -> str:
        """Returns the URL of a file.

        With a ``url_base`` (i.e. "/media/") the URL is the base and the quoted
        name, the same as Django's FileSystemStorage. Use this when the files
        are served by your own site, not by S3. Without a ``url_base``,
        the URL is the presigned URL (see ``presigned_url()``).
        """
        if self.url_base is None:
            return self.presigned_url(name)
        self._key(name)  # Rejects names outside the storage location
        return urljoin(self.url_base, filepath_to_uri(name).lstrip("/"))

    def presigned_url(self, name: str) -> str:
        """Returns a presigned GET URL that expires after ``url_expiry`` seconds."""
        return self.client.generate_presigned_url(
            "get_object",
            Params={"Bucket": self.bucket_name, "Key": self._key(name)},
            ExpiresIn=self.url_expiry,
        )
