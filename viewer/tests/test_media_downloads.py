# pylint: disable=redefined-outer-name,unused-argument,protected-access
"""Tests for structure downloads through the media storage seam (issue 1038).

The download builder (``DownloadStructures``) assembles the archive in a
scratch directory. It used to symlink each media file there, and read the
MOL/SDF files to set their molecule names. Now it queues them and gets them
all with one parallel ``media_storage.fetch()`` (``fetch_queued()``), because
each S3 request is slow. A file "exists" when it is in one listing of the
target's directory. The finished archive is stored as
``downloads/<task_id>/<file>`` and served through nginx (``media_response``).

Each test runs with both media storages (the ``media_backend`` fixture).
"""
import shutil
import tarfile
from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest
from rdkit import Chem

from viewer import download_structures, media_storage
from viewer.download_structures import ArchiveFile, DownloadStructures
from viewer.models import (
    ComputedSet,
    ComputedSetSiteObservation,
    DownloadLinks,
    ExperimentUpload,
    SiteObservation,
)

_TARGET_DIR = "target_loader_data/Mpro_lb1"
_MOL_NO_NAME = "\n" + Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")).split("\n", 1)[1]


@pytest.fixture
def target(user, make_project, make_target):
    target = make_target(make_project("lb1", members=[user]), title="Mpro")
    target.zip_archive = "Mpro_lb1"
    target.save()
    return target


@pytest.fixture
def downloader(media_backend, target, tmp_path):
    tempdir = tmp_path / "download"
    tempdir.mkdir()
    return DownloadStructures(
        tempdir=str(tempdir),
        task=SimpleNamespace(
            request=SimpleNamespace(id="task-1"), update_state=lambda **_: None
        ),
        target=target,
        use_zip=False,
        target_access_string="lb1",
    )


def _archived(downloader) -> dict[str, bytes]:
    root = downloader.temp_path
    return {
        p.relative_to(root).as_posix(): p.read_bytes()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


# DownloadLinks ----------------------------------------------------------------------


def test_get_file_url_is_the_storage_name():
    link = DownloadLinks(task_id="t1", file_url="Mpro.tar.gz")

    assert link.get_file_url() == "downloads/t1/Mpro.tar.gz"
    assert DownloadLinks(task_id="t1").get_file_url() is None
    assert DownloadLinks(file_url="Mpro.tar.gz").get_file_url() is None


# Collecting files ---------------------------------------------------------------------


def test_files_are_copied_into_the_archive(downloader):
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/x0001.pdb", b"ATOM\n")
    archive_file = ArchiveFile(
        name=f"{_TARGET_DIR}/upload_1/x0001.pdb",
        archive_path="aligned_files/x0001/x0001.pdb",
    )

    assert downloader._add_file_to_zip_aligned("x0001", archive_file)
    downloader.fetch_queued()

    assert _archived(downloader) == {"aligned_files/x0001/x0001.pdb": b"ATOM\n"}


def test_mol_files_get_the_molecule_name(downloader):
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/x0001.mol", _MOL_NO_NAME)
    archive_file = ArchiveFile(
        name=f"{_TARGET_DIR}/upload_1/x0001.mol",
        archive_path="aligned_files/x0001/x0001.mol",
    )

    assert downloader._add_file_to_zip_aligned("x0001", archive_file)
    downloader.fetch_queued()

    content = _archived(downloader)["aligned_files/x0001/x0001.mol"].decode()
    assert content.startswith("x0001\n")
    assert content.endswith("$$$$\n\n")


@pytest.mark.parametrize("path", [f"{_TARGET_DIR}/upload_1/missing.pdb", "None", ""])
def test_a_missing_file_gives_an_empty_marker(downloader, path):
    archive_file = ArchiveFile(name=path, archive_path="aligned_files/x0001/x.pdb")

    assert not downloader._add_file_to_zip_aligned("x0001", archive_file)
    downloader.fetch_queued()

    assert _archived(downloader) == {
        "aligned_files/x0001/x.pdb_FILE_NOT_IN_UPLOAD": b""
    }


def test_the_combined_sdf_keeps_the_order(downloader, settings):
    for code in ("b", "a"):
        media_storage.store_bytes(
            f"{_TARGET_DIR}/upload_1/{code}.sdf", f"{code}\n$$$$\n"
        )
    combined = str(downloader.combined_sdf_path)

    for code in ("b", "a"):
        # The builder gives legacy absolute SDF paths (MEDIA_ROOT joined)
        path = str(Path(settings.MEDIA_ROOT) / f"{_TARGET_DIR}/upload_1/{code}.sdf")
        archive_file = ArchiveFile(name=path, archive_path=f"{code}.sdf")
        assert downloader._add_file_to_sdf(combined, archive_file)
    downloader.fetch_queued()

    assert Path(combined).read_text() == "b\n$$$$\na\n$$$$\n"


def _upload(target, user, **fields):
    return ExperimentUpload.objects.create(
        project=target.project,
        target=target,
        file="Mpro.tgz",
        upload_data_dir="upload_1",
        commit_datetime=datetime(2026, 1, 1, tzinfo=timezone.utc),
        committer=user,
        **fields,
    )


def test_yaml_files_are_added(downloader, target, user):
    _upload(target, user, assembly_transforms=f"{_TARGET_DIR}/upload_1/assemblies.yaml")
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/meta_aligner.yaml", b"meta")
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/assemblies.yaml", b"transform")
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/sub/other.yaml", b"not top")

    downloader._yaml_files_zip(target)
    downloader.fetch_queued()

    assert _archived(downloader) == {"yaml_files/upload_1/meta_aligner.yaml": b"meta"}


@pytest.mark.parametrize("soakdb_files", [True, False])
def test_extra_files_are_added(downloader, target, user, soakdb_files):
    _upload(target, user)
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/extra_files/a/notes.txt", b"n")
    media_storage.store_bytes(
        f"{_TARGET_DIR}/upload_1/extra_files/soakdb_1.sqlite", b"s"
    )

    downloader._extra_files_zip(target, soakdb_files=soakdb_files)
    downloader.fetch_queued()

    expected = {"extra_files_1/a/notes.txt": b"n"}
    if soakdb_files:
        expected["extra_files_1/soakdb_1.sqlite"] = b"s"
    assert _archived(downloader) == expected


def test_trans_matrix_files_are_added(downloader, target, user):
    _upload(
        target,
        user,
        neighbourhood_transforms=f"{_TARGET_DIR}/upload_1/neighbourhoods.yaml",
        assembly_transforms=f"{_TARGET_DIR}/upload_1/assemblies.yaml",
    )
    media_storage.store_bytes(f"{_TARGET_DIR}/upload_1/neighbourhoods.yaml", b"n")

    downloader._trans_matrix_files_zip(target)
    downloader.fetch_queued()

    # The transformation matrices go to the root of the archive
    assert _archived(downloader) == {
        "neighbourhoods.yaml": b"n",
        "assemblies.yaml_FILE_NOT_IN_UPLOAD": b"",
    }


def test_compound_set_mols_are_read_from_storage(downloader, target, user):
    mol = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
    media_storage.store_bytes("computed_set_data/good.mol", mol)
    cset = ComputedSet.objects.create(
        name="cset", target=target, owner_user=user, submitted_sdf="hits.sdf"
    )
    for name, mol_file in (("good-hit", "good.mol"), ("missing-hit", "gone.mol")):
        obs = SiteObservation.objects.create(
            smiles="CCO",
            virtual_name=name,
            virtual_ligand_mol=f"computed_set_data/{mol_file}",
        )
        ComputedSetSiteObservation.objects.create(
            computed_set=cset, site_observation=obs
        )

    downloader._compound_sets_zip(target)

    written = _archived(downloader)["virtual_hits/hits.sdf"].decode()
    assert "good-hit" in written
    assert "missing-hit" not in written


# The archive -----------------------------------------------------------------------------


def test_the_archive_is_stored_under_downloads(downloader, tmp_path):
    downloader.write_file("readme", "README.md")

    downloader.compress_and_store("downloads/task-1/Mpro.tar.gz")

    local = tmp_path / "Mpro.tar.gz"
    media_storage.fetch([("downloads/task-1/Mpro.tar.gz", local)])
    with tarfile.open(local) as archive:
        assert "./README.md" in archive.getnames()
    assert not list(downloader.temp_path.glob("*.tar.gz"))


@pytest.mark.skipif(shutil.which("7z") is None, reason="7z is not installed")
def test_a_zip_archive_is_stored_under_downloads(downloader):
    downloader.use_zip = True
    downloader.write_file("readme", "README.md")

    downloader.compress_and_store("downloads/task-1/Mpro.zip")

    assert media_storage.size("downloads/task-1/Mpro.zip") > 0


# Hard erase ---------------------------------------------------------------------------------


@pytest.mark.django_db
def test_hard_erase_deletes_the_download_directory(media_backend):
    long_ago = datetime.now(timezone.utc) - timedelta(days=30)
    DownloadLinks.objects.create(
        task_id="dead",
        file_url="Mpro.tar.gz",
        create_date=long_ago,
        expired_date=long_ago,
    )
    media_storage.store_bytes("downloads/dead/Mpro.tar.gz", b"x")
    media_storage.store_bytes("downloads/dead_2/Mpro.tar.gz", b"keep")

    download_structures.hard_erase_out_of_date_download_records()

    assert [o.name for o in media_storage.list_prefix("downloads")] == [
        "downloads/dead_2/Mpro.tar.gz"
    ]
    assert DownloadLinks.objects.get(task_id="dead").deleted


# Serving ------------------------------------------------------------------------------------------


@pytest.fixture
def stored_download(media_backend, target, user):
    DownloadLinks.objects.create(
        task_id="t1",
        file_url="Mpro.tar.gz",
        target=target,
        user=user,
        create_date=datetime.now(timezone.utc),
    )
    media_storage.store_bytes("downloads/t1/Mpro.tar.gz", b"archive")


@pytest.mark.parametrize(
    "file_url", ["downloads/t1/Mpro.tar.gz", "/code/media/downloads/t1/Mpro.tar.gz"]
)
def test_a_download_is_served_through_nginx(
    stored_download, media_backend, authenticated_client, file_url
):
    response = authenticated_client.get(
        f"/api/download_structures/?file_url={file_url}"
    )

    assert response.status_code == 200
    location = "/media-s3/" if media_backend == "s3" else "/media-local/"
    assert response["X-Accel-Redirect"] == location + "downloads/t1/Mpro.tar.gz"
    assert response["Content-Disposition"] == 'attachment; filename="Mpro.tar.gz"'


def test_head_gives_the_download_size(stored_download, authenticated_client):
    response = authenticated_client.head(
        "/api/download_structures/?file_url=downloads/t1/Mpro.tar.gz"
    )

    assert response.status_code == 200
    assert response["Content-Length"] == str(len(b"archive"))
    assert "X-Accel-Redirect" not in response


def test_an_unknown_download_is_not_found(stored_download, authenticated_client):
    response = authenticated_client.get(
        "/api/download_structures/?file_url=downloads/other/Mpro.tar.gz"
    )

    assert response.status_code == 404
