"""Access to the backend's media files, for either media storage (issue 1038).

MEDIA_STORAGE selects the storage: the local filesystem (under MEDIA_ROOT) or
an S3 bucket (see ``frag_s3_storage``). All code that reads, writes, moves,
lists, deletes or serves media files uses this module, so that it does not need
to know which storage is in use. Do not use MEDIA_ROOT, a FieldFile's ``.path``
or ``default_storage`` anywhere else (``viewer/tests/test_media_seam.py`` checks).

Names are storage-relative POSIX names (i.e. ``target_loader_data/Mpro/a.pdb``),
the same as a FileField's ``name``. Functions that transfer many files use
MEDIA_TRANSFER_THREADS threads, because one S3 request takes 30-75 ms.
Work that needs local files (tar, RDKit, 7z) uses scratch copies: see
``fetch()``, ``local_copy()`` and ``store_tree()``. Celery tasks must pass
storage names to each other, never local paths, because each task can run on
a different worker.

Media files are served by nginx, not by Django. A view checks that the user
can see the file, then returns ``media_response()``. That response has an
``X-Accel-Redirect`` to an internal nginx location (see "django_nginx.conf"):

- ``/media-local/<name>``: nginx reads the file from MEDIA_ROOT.
- ``/media-s3/<name>``: nginx proxies the presigned S3 URL in the
  ``X-Media-Url`` header. The browser never sees that URL, so the file stays
  on the same origin (the page is cross-origin isolated) and the URL
  cannot be shared.
"""

import os
import posixpath
import shutil
import tempfile
import uuid
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, NamedTuple, TypeVar
from urllib.parse import quote

from django.conf import settings
from django.core.exceptions import SuspiciousFileOperation
from django.core.files.base import ContentFile
from django.core.files.storage import FileSystemStorage, Storage, storages
from django.core.files.uploadedfile import UploadedFile
from django.http import Http404, HttpRequest, HttpResponse
from django.utils.text import get_valid_filename

from frag_s3_storage import S3Storage

# The internal nginx locations, see "django_nginx.conf".
_LOCAL_LOCATION = "/media-local/"
_S3_LOCATION = "/media-s3/"

_T = TypeVar("_T")
_R = TypeVar("_R")


class StoredObject(NamedTuple):
    """A file in the media storage."""

    name: str
    size: int
    modified: datetime


def _storage() -> Storage:
    # Not cached: tests change the STORAGES setting.
    return storages["default"]


def _local_path(name: str) -> Path:
    """The local path of a file in FILESYSTEM media storage."""
    storage = _storage()
    assert isinstance(storage, FileSystemStorage)
    return Path(storage.path(name))


def _parallel(func: Callable[[_T], _R], items: Iterable[_T]) -> list[_R]:
    """Calls func for each item, in MEDIA_TRANSFER_THREADS threads.

    The results are in the order of the items. The first exception is raised.
    """
    storage = _storage()
    if isinstance(storage, S3Storage):
        # Make the (thread-safe) client now, not in each thread.
        _ = storage.client
    threads = max(1, settings.MEDIA_TRANSFER_THREADS)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        return list(executor.map(func, items))


# Names -----------------------------------------------------------------------


def is_s3() -> bool:
    """True if the media files are in S3."""
    return isinstance(_storage(), S3Storage)


def to_name(value: Any) -> str:
    """Returns the clean storage name of a value (a name, Path or FieldFile).

    A legacy absolute path under MEDIA_ROOT (i.e. "/code/media/a/b.sdf")
    becomes a storage name ("a/b.sdf"). A name that is outside the
    media storage raises SuspiciousFileOperation.
    """
    text = str(value)
    media_root = str(settings.MEDIA_ROOT).rstrip("/") + "/"
    if text.startswith(media_root):
        text = text[len(media_root) :]
    clean = posixpath.normpath(text) if text else ""
    if clean in ("", ".", "..") or clean.startswith(("/", "../")) or "\\" in clean:
        raise SuspiciousFileOperation(f"'{value}' is not a media file name")
    return clean


# Read and write --------------------------------------------------------------


def exists(name: str) -> bool:
    return _storage().exists(to_name(name))


def size(name: str) -> int:
    """The size of the file in bytes. Raises FileNotFoundError for a missing file."""
    return _storage().size(to_name(name))


def read_bytes(name: str) -> bytes:
    """The content of the file. Raises FileNotFoundError for a missing file."""
    with _storage().open(to_name(name), "rb") as file:
        return file.read()


def read_text(name: str, encoding: str = "utf-8") -> str:
    """The content of the file as text. Raises FileNotFoundError for a missing file."""
    return read_bytes(name).decode(encoding)


def store_bytes(name: str, data: bytes | str) -> str:
    """Stores the data (text as UTF-8) with exactly this name.

    An existing file with the name is replaced. Returns the name.
    """
    name = to_name(name)
    if isinstance(data, str):
        data = data.encode("utf-8")
    storage = _storage()
    # Storage.save() adds a suffix to a name that exists, so delete it first.
    storage.delete(name)
    saved = storage.save(name, ContentFile(data))
    if saved != name:
        raise RuntimeError(f"'{name}' was stored as '{saved}'")
    return name


def store_file(local_path: str | os.PathLike, name: str) -> str:
    """Stores a local file with exactly this name, and replaces any existing file.

    Returns the name.
    """
    name = to_name(name)
    storage = _storage()
    if isinstance(storage, S3Storage):
        storage.upload_file(local_path, name)
    else:
        target = _local_path(name)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(local_path, target)
    return name


def store_upload(uploaded: UploadedFile, prefix: str) -> str:
    """Stores an uploaded file as ``<prefix>/<uuid>/<safe file name>``.

    Each upload gets a new directory, so two uploads of the same file do not
    replace each other. Returns the name.
    """
    file_name = get_valid_filename(PurePosixPath(str(uploaded.name)).name)
    name = to_name(f"{prefix}/{uuid.uuid4().hex}/{file_name}")
    if hasattr(uploaded, "temporary_file_path"):
        return store_file(uploaded.temporary_file_path(), name)
    uploaded.seek(0)
    saved = _storage().save(name, uploaded)
    if saved != name:
        raise RuntimeError(f"'{name}' was stored as '{saved}'")
    return name


# Trees and scratch copies ---------------------------------------------------


def store_tree(local_dir: str | os.PathLike, prefix: str) -> int:
    """Stores every file below a local directory under the prefix, in parallel.

    The names are the prefix and each file's path relative to the directory.
    Existing files are replaced. Returns the number of files.
    """
    root = Path(local_dir)
    prefix = to_name(prefix)
    pairs = [
        (path, f"{prefix}/{path.relative_to(root).as_posix()}")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    ]
    _parallel(lambda pair: store_file(*pair), pairs)
    return len(pairs)


def fetch(pairs: Iterable[tuple[str, str | os.PathLike]]) -> list[str]:
    """Copies files to local paths, in parallel. Each pair is (name, local path).

    Missing directories are made. Returns the names of the missing files,
    in place of an exception, so that the caller decides what to do.
    """
    storage = _storage()

    def _fetch_one(pair: tuple[str, str | os.PathLike]) -> str | None:
        name, local_path = to_name(pair[0]), Path(pair[1])
        local_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            if isinstance(storage, S3Storage):
                storage.download_file(name, local_path)
            else:
                shutil.copyfile(_local_path(name), local_path)
        except FileNotFoundError:
            return name
        return None

    return [name for name in _parallel(_fetch_one, list(pairs)) if name]


@contextmanager
def local_copy(name: str) -> Iterator[Path]:
    """Gives a scratch copy of the file, with the same file name.

    The copy (and its temporary directory) is removed at the end. Changes to
    the copy do not change the stored file. Raises FileNotFoundError for a
    missing file.
    """
    with tempfile.TemporaryDirectory() as scratch:
        path = Path(scratch) / PurePosixPath(to_name(name)).name
        if fetch([(name, path)]):
            raise FileNotFoundError(name)
        yield path


# Move and delete ------------------------------------------------------------


def move(source: str, target: str) -> None:
    """Renames a file. Raises FileNotFoundError for a missing source file."""
    source, target = to_name(source), to_name(target)
    storage = _storage()
    if isinstance(storage, S3Storage):
        # S3 has no rename.
        storage.copy(source, target)
        storage.delete(source)
    else:
        target_path = _local_path(target)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        os.replace(_local_path(source), target_path)


def delete(name: str) -> None:
    """Deletes the file. A missing file is not an error."""
    _storage().delete(to_name(name))


def delete_many(names: Iterable[str]) -> int:
    """Deletes the files (S3: 1000 for each request). Missing files are not errors.

    Returns the number of names. Raises OSError if a file is not deleted.
    """
    clean = [to_name(name) for name in names]
    storage = _storage()
    if isinstance(storage, S3Storage):
        failed = storage.delete_many(clean)
        if failed:
            raise OSError(
                f"{len(failed)} media files are not deleted, i.e. {failed[0]}"
            )
    else:
        for name in clean:
            storage.delete(name)
    return len(clean)


def delete_prefix(prefix: str) -> int:
    """Deletes every file below a directory. Returns the number of files.

    ``"t/Mpro"`` deletes ``t/Mpro/...`` but not ``t/Mpro_2/...``.
    The whole storage (``""``, ``"."`` or ``"/"``) is refused (ValueError).
    """
    if prefix.strip() in ("", ".", "/"):
        raise ValueError("delete_prefix() does not delete the whole media storage")
    prefix = to_name(prefix)
    names = [stored.name for stored in list_prefix(prefix)]
    delete_many(names)
    if not is_s3():
        # Also remove the (now empty) directories.
        shutil.rmtree(_local_path(prefix), ignore_errors=True)
    return len(names)


# Listing ---------------------------------------------------------------------


def list_prefix(prefix: str) -> Iterator[StoredObject]:
    """Yields every file below a directory, at any depth.

    ``"t/Mpro"`` lists ``t/Mpro/...`` but not ``t/Mpro_2/...``.
    ``""`` lists the whole storage. A missing directory gives nothing.
    """
    prefix = to_name(prefix) if prefix else ""
    storage = _storage()
    if isinstance(storage, S3Storage):
        for name, file_size, modified in storage.iter_objects(
            f"{prefix}/" if prefix else ""
        ):
            yield StoredObject(name, file_size, modified)
        return

    assert isinstance(storage, FileSystemStorage)
    root = Path(storage.location)
    base = root / prefix if prefix else root
    if not base.is_dir():
        return
    for directory, _, files in os.walk(base):
        for file_name in sorted(files):
            path = Path(directory) / file_name
            try:
                stat = path.stat()
            except FileNotFoundError:
                continue  # A broken symbolic link
            yield StoredObject(
                path.relative_to(root).as_posix(),
                stat.st_size,
                datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc),
            )


def list_dirs(prefix: str) -> list[str]:
    """The names of the directories in a directory (not their full names), sorted.

    A missing directory gives an empty list.
    """
    prefix = to_name(prefix) if prefix else ""
    try:
        directories, _ = _storage().listdir(prefix)
    except FileNotFoundError:
        return []
    return sorted(directories)


# Serving ---------------------------------------------------------------------


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
    try:
        name = to_name(name)
    except SuspiciousFileOperation as exc:
        raise Http404(f"'{name}' is not a media file") from exc
    storage = _storage()

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
