"""Tests for the orphaned-download recovery (``viewer.download_restore``).

Migration 0153 wiped every ``DownloadLinks`` row on the deploy that applied it,
leaving the archives under ``media/downloads/`` with nothing pointing at them.
The recovery reads the original POST body back out of the README.md inside each
archive (written by ``DownloadStructures._build_readme()``) and rebuilds the row.

The archives here are built the way the real ones are - from a directory, with
the README at the root - in both flavours (``.tar.gz`` and ``.zip``), so the
member-name handling is exercised for real rather than mocked.
"""
# Fixtures legitimately reuse their names as arguments, and `db` is requested
# for its side effect only - both are standard pytest patterns pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument

import json
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pytest

from viewer.download_restore import (
    DOWNLOADS_SUBDIR,
    read_original_search,
    restore,
    scan_downloads,
)
from viewer.download_structures import KEEP_UNTIL_DURATION
from viewer.models import DownloadLinks

TASK_ID = "4c3afc69-bca9-4fb1-a76e-56c85a85899f"


def _search(target="Mpro", tas="lb32627", static_link=True) -> dict:
    """A POST body of the shape the frontend sends."""
    return {
        "target_name": target,
        "target_access_string": tas,
        "proteins": "Mpro-x0104_A,Mpro-x0107_A",
        "all_aligned_structures": True,
        "metadata_info": True,
        "static_link": static_link,
    }


def _readme(original_search: dict) -> str:
    """The README as _build_readme() writes it - the JSON on its own fenced line."""
    return (
        "# Documentation for the downloaded zipfile\n"
        "## Download details\n"
        "\n### Download command (JSON)\n"
        "JSON command sent from front-end to backend to generate the download.\n\n"
        f"```{json.dumps(original_search)}```\n\n"
        "## Files included\n- README.md\n"
    )


@pytest.fixture
def now_utc() -> datetime:
    return datetime.now(timezone.utc)


@pytest.fixture
def downloads(settings, tmp_path) -> Path:
    """Point MEDIA_ROOT at a temp dir and return its downloads/ subdirectory."""
    media_root = tmp_path / "media"
    downloads_dir = media_root / DOWNLOADS_SUBDIR
    downloads_dir.mkdir(parents=True)
    settings.MEDIA_ROOT = str(media_root)
    return downloads_dir


def _make_archive(
    downloads: Path,
    *,
    task_id: str = TASK_ID,
    name: str = "Mpro.tar.gz",
    original_search: dict | None = None,
    with_readme: bool = True,
) -> Path:
    """Lay down media/downloads/<task_id>/<name> the way create_download() does."""
    directory = downloads / task_id
    directory.mkdir(parents=True, exist_ok=True)

    content = directory / "content"
    content.mkdir(exist_ok=True)
    (content / "metadata.csv").write_text("code,smiles\n")
    if with_readme:
        search = _search() if original_search is None else original_search
        (content / "README.md").write_text(_readme(search))

    archive = directory / name
    if name.endswith(".zip"):
        with zipfile.ZipFile(archive, "w") as zip_file:
            for path in sorted(content.rglob("*")):
                zip_file.write(path, path.relative_to(content).as_posix())
    else:
        with tarfile.open(archive, "w:gz") as tar_file:
            # tar is invoked as `tar -C <dir> ... .`, so members are "./<name>".
            tar_file.add(content, arcname=".")
    return archive


# --------------------------------------------------------------------------- #
# Reading the search back out of an archive
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["Mpro.tar.gz", "Mpro.zip"])
def test_original_search_is_read_from_either_archive_flavour(downloads, name):
    archive = _make_archive(downloads, name=name)

    assert read_original_search(archive) == _search()


def test_archive_without_readme_yields_no_search(downloads):
    archive = _make_archive(downloads, with_readme=False)

    assert read_original_search(archive) is None


def test_unreadable_archive_yields_no_search(downloads):
    directory = downloads / TASK_ID
    directory.mkdir(parents=True)
    archive = directory / "Mpro.tar.gz"
    archive.write_text("not a tarball")

    assert read_original_search(archive) is None


# --------------------------------------------------------------------------- #
# Scanning
# --------------------------------------------------------------------------- #


@pytest.mark.django_db
def test_scan_finds_orphan_and_resolves_its_target(
    downloads, make_project, make_target
):
    project = make_project("lb32627")
    target = make_target(project, title="Mpro")
    _make_archive(downloads)

    report = scan_downloads()

    assert len(report.orphans) == 1
    orphan = report.orphans[0]
    assert orphan.task_id == TASK_ID
    assert orphan.file_url == "Mpro.tar.gz"
    assert orphan.static_link is True
    assert orphan.target == target
    assert report.static_orphans == [orphan]


@pytest.mark.django_db
def test_scan_ignores_directories_that_still_have_a_record(downloads, now_utc):
    _make_archive(downloads)
    DownloadLinks.objects.create(
        task_id=TASK_ID, file_url="Mpro.tar.gz", create_date=now_utc
    )

    report = scan_downloads()

    assert report.orphans == []
    assert report.referenced == 1


@pytest.mark.django_db
def test_scan_skips_directories_without_an_archive(downloads):
    (downloads / "empty-dir").mkdir()

    report = scan_downloads()

    assert report.orphans == []
    assert [reason for _, reason in report.skipped] == ["no .zip/.tar.gz archive"]


@pytest.mark.django_db
def test_scan_separates_dynamic_and_unparsed_orphans(downloads):
    _make_archive(downloads, task_id="static-1", original_search=_search())
    _make_archive(
        downloads, task_id="dynamic-1", original_search=_search(static_link=False)
    )
    _make_archive(downloads, task_id="unparsed-1", with_readme=False)

    report = scan_downloads()

    assert [o.task_id for o in report.static_orphans] == ["static-1"]
    assert [o.task_id for o in report.dynamic_orphans] == ["dynamic-1"]
    assert [o.task_id for o in report.unparsed_orphans] == ["unparsed-1"]


# --------------------------------------------------------------------------- #
# Restoring
# --------------------------------------------------------------------------- #


@pytest.mark.django_db
def test_restore_rebuilds_the_record_the_lookup_needs(
    downloads, make_project, make_target
):
    """The restored row is exactly what DownloadStructuresView.create() looks up:
    task_id + file_url, with static_link set."""
    project = make_project("lb32627")
    target = make_target(project, title="Mpro")
    archive = _make_archive(downloads)

    result = restore(scan_downloads().orphans)

    assert len(result.restored) == 1
    link = DownloadLinks.objects.get(task_id=TASK_ID)
    assert link.file_url == "Mpro.tar.gz"
    assert link.static_link is True
    assert link.target == target
    assert link.get_file_url() == str(archive)
    assert link.original_search == _search()
    # Derived through the view's own serializer, so the flags match.
    assert link.protein_params["bound_file"] is True
    assert link.other_params["metadata_info"] is True
    # create_date falls back to the archive's mtime.
    assert link.create_date.astimezone(timezone.utc).timestamp() == pytest.approx(
        archive.stat().st_mtime, abs=1
    )


@pytest.mark.django_db
def test_restored_static_link_is_left_alone_by_housekeeping(downloads):
    """No keep_zip_until and static_link=True: all three cleanup jobs skip it."""
    _make_archive(downloads)

    restore(scan_downloads().orphans)

    link = DownloadLinks.objects.get(task_id=TASK_ID)
    assert link.keep_zip_until is None
    assert link.expired_date is None


@pytest.mark.django_db
def test_proteins_is_left_null_so_the_row_is_never_reused(downloads):
    """proteins holds SiteObservation pks resolved at request time - meaningless
    after a re-upload - so a restored row must never satisfy the reuse query."""
    _make_archive(downloads)

    restore(scan_downloads().orphans)

    assert DownloadLinks.objects.get(task_id=TASK_ID).proteins is None


@pytest.mark.django_db
def test_dynamic_links_are_skipped_unless_asked_for(downloads):
    _make_archive(downloads, original_search=_search(static_link=False))

    assert restore(scan_downloads().orphans).restored == []
    assert not DownloadLinks.objects.exists()


@pytest.mark.django_db
def test_restored_dynamic_link_keeps_its_original_retention(downloads):
    archive = _make_archive(downloads, original_search=_search(static_link=False))

    result = restore(scan_downloads().orphans, include_dynamic=True)

    assert len(result.restored) == 1
    link = DownloadLinks.objects.get(task_id=TASK_ID)
    assert link.static_link is False
    expected = result.restored[0].create_date + KEEP_UNTIL_DURATION
    assert link.keep_zip_until == expected
    assert archive.is_file()


@pytest.mark.django_db
def test_unparsed_archives_are_skipped_unless_asked_for(downloads):
    _make_archive(downloads, with_readme=False)

    assert restore(scan_downloads().orphans).restored == []

    result = restore(scan_downloads().orphans, unparsed_as_static=True)

    assert len(result.restored) == 1
    link = DownloadLinks.objects.get(task_id=TASK_ID)
    assert link.static_link is True
    assert link.original_search is None
    assert link.file_url == "Mpro.tar.gz"


@pytest.mark.django_db
def test_restore_is_idempotent(downloads):
    _make_archive(downloads)

    restore(scan_downloads().orphans)
    restore(scan_downloads().orphans)

    assert DownloadLinks.objects.filter(task_id=TASK_ID).count() == 1


@pytest.mark.django_db
def test_restore_survives_a_record_appearing_mid_run(downloads, now_utc):
    """A concurrent download claiming the task_id is reported, not raised."""
    _make_archive(downloads)
    orphans = scan_downloads().orphans
    DownloadLinks.objects.create(
        task_id=TASK_ID, file_url="Mpro.tar.gz", create_date=now_utc
    )

    result = restore(orphans)

    assert result.restored == []
    assert [message for _, message in result.failed] == [
        "a record now exists for this task"
    ]
    assert DownloadLinks.objects.filter(task_id=TASK_ID).count() == 1
