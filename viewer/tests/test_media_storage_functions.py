# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the media storage functions in ``viewer.media_storage`` (issue 1038).

Code that reads, writes, moves, lists or deletes media files uses these
functions, never MEDIA_ROOT or a FieldFile's ``.path``. Each test runs with
both media storages (the ``media_backend`` fixture): the filesystem and S3.
Names are storage-relative POSIX names (i.e. ``target_loader_data/Mpro/a.pdb``).
"""
import datetime
import re
from pathlib import Path

import pytest
from django.core.exceptions import SuspiciousFileOperation
from django.core.files.uploadedfile import SimpleUploadedFile, TemporaryUploadedFile

from viewer import media_storage

_NAME = "target_loader_data/Mpro/aligned/a b#1.pdb"


# Names -----------------------------------------------------------------------


def test_is_s3_matches_the_media_storage(media_backend):
    assert media_storage.is_s3() == (media_backend == "s3")


@pytest.mark.parametrize(
    "value, name",
    [
        ("pdbs/a.pdb", "pdbs/a.pdb"),
        ("pdbs//sub/./a.pdb", "pdbs/sub/a.pdb"),
        (Path("pdbs/a.pdb"), "pdbs/a.pdb"),
    ],
)
def test_to_name_returns_a_clean_storage_name(media_backend, value, name):
    assert media_storage.to_name(value) == name


def test_to_name_strips_a_legacy_media_root_prefix(media_backend, settings):
    value = str(Path(settings.MEDIA_ROOT) / "computed_set_data" / "a.sdf")

    assert media_storage.to_name(value) == "computed_set_data/a.sdf"


@pytest.mark.parametrize(
    "value", ["", ".", "..", "../a", "pdbs/../../a", "/etc/passwd"]
)
def test_to_name_rejects_names_outside_the_media_storage(media_backend, value):
    with pytest.raises(SuspiciousFileOperation):
        media_storage.to_name(value)


# Read and write --------------------------------------------------------------


def test_store_bytes_then_read_returns_the_same_content(media_backend):
    assert media_storage.store_bytes(_NAME, b"ATOM\n") == _NAME

    assert media_storage.read_bytes(_NAME) == b"ATOM\n"
    assert media_storage.read_text(_NAME) == "ATOM\n"
    assert media_storage.exists(_NAME)
    assert media_storage.size(_NAME) == 5


def test_store_bytes_stores_text_as_utf8(media_backend):
    media_storage.store_bytes("a.txt", "µ")

    assert media_storage.read_bytes("a.txt") == "µ".encode("utf-8")


def test_store_bytes_replaces_an_existing_file(media_backend):
    media_storage.store_bytes("a.txt", b"old")

    assert media_storage.store_bytes("a.txt", b"new") == "a.txt"

    assert media_storage.read_bytes("a.txt") == b"new"
    assert media_storage.list_dirs("") == []
    assert [o.name for o in media_storage.list_prefix("")] == ["a.txt"]


def test_a_missing_file(media_backend):
    assert not media_storage.exists("missing.txt")
    with pytest.raises(FileNotFoundError):
        media_storage.read_bytes("missing.txt")
    with pytest.raises(FileNotFoundError):
        media_storage.read_text("missing.txt")
    with pytest.raises(FileNotFoundError):
        media_storage.size("missing.txt")


def test_store_file_stores_a_local_file(media_backend, tmp_path):
    local = tmp_path / "local.pdb"
    local.write_bytes(b"ATOM\n")

    assert media_storage.store_file(local, _NAME) == _NAME

    assert media_storage.read_bytes(_NAME) == b"ATOM\n"


def test_store_file_replaces_an_existing_file(media_backend, tmp_path):
    media_storage.store_bytes("a.txt", b"old")
    local = tmp_path / "a.txt"
    local.write_bytes(b"new")

    media_storage.store_file(local, "a.txt")

    assert media_storage.read_bytes("a.txt") == b"new"


_UPLOAD_NAME = re.compile(r"^tmp/uploads/[0-9a-f]{32}/Mpro_bundle.tar.gz$")


def test_store_upload_stores_an_in_memory_upload_under_a_new_directory(media_backend):
    upload = SimpleUploadedFile("Mpro bundle.tar.gz", b"bundle")

    name = media_storage.store_upload(upload, "tmp/uploads")

    assert _UPLOAD_NAME.match(name)
    assert media_storage.read_bytes(name) == b"bundle"


def test_store_upload_stores_a_temporary_file_upload(media_backend):
    upload = TemporaryUploadedFile("Mpro bundle.tar.gz", "application/gzip", 6, None)
    upload.write(b"bundle")
    upload.seek(0)

    name = media_storage.store_upload(upload, "tmp/uploads")

    assert _UPLOAD_NAME.match(name)
    assert media_storage.read_bytes(name) == b"bundle"


def test_store_upload_gives_two_uploads_different_names(media_backend):
    first = media_storage.store_upload(SimpleUploadedFile("a.txt", b"1"), "tmp")
    second = media_storage.store_upload(SimpleUploadedFile("a.txt", b"2"), "tmp")

    assert first != second


# Trees and scratch copies ---------------------------------------------------


def _make_tree(root: Path) -> None:
    (root / "aligned" / "x0001").mkdir(parents=True)
    (root / "aligned" / "x0001" / "a b#1.pdb").write_bytes(b"pdb")
    (root / "aligned" / "x0001" / "a.map").write_bytes(b"map" * 100)
    (root / "meta.csv").write_bytes(b"csv")


def test_store_tree_stores_every_file_under_the_prefix(media_backend, tmp_path):
    local = tmp_path / "local"
    _make_tree(local)

    count = media_storage.store_tree(local, "target_loader_data/Mpro")

    assert count == 3
    assert sorted(o.name for o in media_storage.list_prefix("target_loader_data")) == [
        "target_loader_data/Mpro/aligned/x0001/a b#1.pdb",
        "target_loader_data/Mpro/aligned/x0001/a.map",
        "target_loader_data/Mpro/meta.csv",
    ]
    assert media_storage.read_bytes("target_loader_data/Mpro/meta.csv") == b"csv"


def test_store_tree_works_with_one_thread(media_backend, tmp_path, settings):
    settings.MEDIA_TRANSFER_THREADS = 1
    local = tmp_path / "local"
    _make_tree(local)

    assert media_storage.store_tree(local, "t") == 3


def test_fetch_copies_files_to_local_paths(media_backend, tmp_path):
    media_storage.store_bytes("a.txt", b"a")
    media_storage.store_bytes(_NAME, b"pdb")
    target = tmp_path / "out"

    missing = media_storage.fetch(
        [("a.txt", target / "a.txt"), (_NAME, target / "deep" / "dir" / "b.pdb")]
    )

    assert missing == []
    assert (target / "a.txt").read_bytes() == b"a"
    assert (target / "deep" / "dir" / "b.pdb").read_bytes() == b"pdb"


def test_fetch_returns_the_missing_names(media_backend, tmp_path):
    media_storage.store_bytes("a.txt", b"a")

    missing = media_storage.fetch(
        [("a.txt", tmp_path / "a.txt"), ("missing.txt", tmp_path / "m.txt")]
    )

    assert missing == ["missing.txt"]
    assert (tmp_path / "a.txt").exists()
    assert not (tmp_path / "m.txt").exists()


def test_local_copy_gives_a_scratch_copy_that_is_removed_after_use(media_backend):
    media_storage.store_bytes(_NAME, b"pdb")

    with media_storage.local_copy(_NAME) as path:
        assert path.name == "a b#1.pdb"
        assert path.read_bytes() == b"pdb"
        path.write_bytes(b"changed")  # A scratch copy, not the stored file

    assert not path.exists()
    assert not path.parent.exists()
    assert media_storage.read_bytes(_NAME) == b"pdb"


def test_local_copy_of_a_missing_file_raises_file_not_found(media_backend):
    with pytest.raises(FileNotFoundError):
        with media_storage.local_copy("missing.txt"):
            pass


# Move and delete ------------------------------------------------------------


def test_move_renames_a_file(media_backend):
    media_storage.store_bytes("tmp/a.txt", b"a")

    media_storage.move("tmp/a.txt", "computed_set_data/b.txt")

    assert not media_storage.exists("tmp/a.txt")
    assert media_storage.read_bytes("computed_set_data/b.txt") == b"a"


def test_move_of_a_missing_file_raises_file_not_found(media_backend):
    with pytest.raises(FileNotFoundError):
        media_storage.move("missing.txt", "b.txt")


def test_delete_removes_a_file_and_ignores_a_missing_file(media_backend):
    media_storage.store_bytes("a.txt", b"a")

    media_storage.delete("a.txt")
    media_storage.delete("a.txt")

    assert not media_storage.exists("a.txt")


def test_delete_many_removes_the_files(media_backend):
    for name in ("a.txt", "dir/b.txt", "dir/c.txt"):
        media_storage.store_bytes(name, b"x")

    media_storage.delete_many(["a.txt", "dir/b.txt", "missing.txt"])

    assert [o.name for o in media_storage.list_prefix("")] == ["dir/c.txt"]


def test_delete_prefix_removes_one_directory_only(media_backend):
    for name in ("t/Mpro/a.txt", "t/Mpro/sub/b.txt", "t/Mpro_2/c.txt", "t/Mpro.txt"):
        media_storage.store_bytes(name, b"x")

    assert media_storage.delete_prefix("t/Mpro") == 2

    assert sorted(o.name for o in media_storage.list_prefix("t")) == [
        "t/Mpro.txt",
        "t/Mpro_2/c.txt",
    ]
    assert media_storage.list_dirs("t") == ["Mpro_2"]


def test_delete_prefix_of_a_missing_directory_deletes_nothing(media_backend):
    assert media_storage.delete_prefix("missing") == 0


@pytest.mark.parametrize("prefix", ["", ".", "/", "../t"])
def test_delete_prefix_refuses_the_whole_storage(media_backend, prefix):
    media_storage.store_bytes("a.txt", b"x")

    with pytest.raises((ValueError, SuspiciousFileOperation)):
        media_storage.delete_prefix(prefix)

    assert media_storage.exists("a.txt")


# Listing ---------------------------------------------------------------------


def test_list_prefix_lists_every_file_below_a_directory(media_backend):
    media_storage.store_bytes("t/Mpro/a.txt", b"abc")
    media_storage.store_bytes("t/Mpro/sub/b.txt", b"x")
    media_storage.store_bytes("t/Mpro_2/c.txt", b"x")

    objects = sorted(media_storage.list_prefix("t/Mpro"))

    assert [(o.name, o.size) for o in objects] == [
        ("t/Mpro/a.txt", 3),
        ("t/Mpro/sub/b.txt", 1),
    ]
    assert all(o.modified.tzinfo is not None for o in objects)
    assert all(isinstance(o.modified, datetime.datetime) for o in objects)


def test_list_prefix_of_a_missing_directory_is_empty(media_backend):
    assert list(media_storage.list_prefix("missing")) == []


def test_list_dirs_lists_the_directories_in_a_directory(media_backend):
    media_storage.store_bytes("t/Mpro/a.txt", b"x")
    media_storage.store_bytes("t/A71EV2A/sub/b.txt", b"x")
    media_storage.store_bytes("t/top.txt", b"x")

    assert media_storage.list_dirs("t") == ["A71EV2A", "Mpro"]


def test_list_dirs_of_a_missing_directory_is_empty(media_backend):
    assert media_storage.list_dirs("missing") == []
