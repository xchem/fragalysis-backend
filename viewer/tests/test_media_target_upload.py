# pylint: disable=redefined-outer-name,unused-argument
"""Tests for target (LHS) uploads through the media storage seam (issue 1038).

The upload view stores the bundle (and an optional curation file) as
``tmp/uploads/<uuid>/<file>`` and gives the Celery task their storage names.
The loader copies the bundle to scratch, extracts it there, and at the end of
its transaction stores the upload's files under
``target_loader_data/<target>/<upload_N>/`` and moves the bundle next to them.

Each storage test runs with both media storages (the ``media_backend`` fixture).
"""
import hashlib
import io
import re
import tarfile
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from django.core.management import call_command
from rdkit import Chem

from frag_s3_storage import S3UploadError
from viewer import media_storage
from viewer.target_loader import (
    TargetLoader,
    UploadReport,
    delete_stored_upload_files,
    load_target,
    read_stored_mol,
    report_store_failure,
    store_upload_files,
)
from viewer.utils import sanitize_directory_name

_MOL_BLOCK = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
_UPLOAD_NAME = re.compile(r"^tmp/uploads/[0-9a-f]{32}/(?P<file>[^/]+)$")


def _upload_file_name(name: str) -> str:
    """The file name in an upload's storage name (tmp/uploads/<uuid>/<file>)."""
    match = _UPLOAD_NAME.match(name)
    assert match, name
    return match["file"]


# The upload view --------------------------------------------------------------


@pytest.fixture
def upload_view(authenticated_client, mocker):
    """The upload endpoint, with authorisation, Celery and the task mocked."""
    mocker.patch("viewer.views._upload_tas_authorisation_response", return_value=None)
    mocker.patch("viewer.views.Celery")
    delay = mocker.patch("viewer.views.task_load_target.delay")
    delay.return_value.task_id = "12345678-1234-1234-1234-123456789abc"

    def _post(**data):
        return authenticated_client.post(
            "/api/upload_target_experiments/",
            {"target_access_string": "lb12345-1", **data},
            format="multipart",
        )

    return _post, delay


def test_upload_stores_the_bundle_and_passes_its_name(media_backend, upload_view):
    post, delay = upload_view

    response = post(file=SimpleUploadedFile("Mpro.tgz", b"bundle"))

    assert response.status_code == 202
    kwargs = delay.call_args.kwargs
    assert _upload_file_name(kwargs["data_bundle"]) == "Mpro.tgz"
    assert media_storage.read_bytes(kwargs["data_bundle"]) == b"bundle"
    assert kwargs["curation_file"] is None


def test_upload_stores_the_curation_file(media_backend, upload_view):
    post, delay = upload_view

    post(
        file=SimpleUploadedFile("Mpro.tgz", b"bundle"),
        curation_file=SimpleUploadedFile("curation.xlsx", b"xlsx"),
    )

    curation = delay.call_args.kwargs["curation_file"]
    assert _upload_file_name(curation) == "curation.xlsx"
    assert media_storage.read_bytes(curation) == b"xlsx"


def test_upload_with_a_good_checksum_is_accepted(media_backend, upload_view):
    post, _ = upload_view

    response = post(
        file=SimpleUploadedFile("Mpro.tgz", b"bundle"),
        sha256checksum=hashlib.sha256(b"bundle").hexdigest(),
    )

    assert response.status_code == 202


def test_upload_with_a_bad_checksum_stores_nothing(media_backend, upload_view):
    post, delay = upload_view

    response = post(file=SimpleUploadedFile("Mpro.tgz", b"bundle"), sha256checksum="0")

    assert response.status_code == 400
    assert delay.call_count == 0
    assert list(media_storage.list_prefix("tmp")) == []


# The load_target command -------------------------------------------------------


def test_load_target_command_stores_the_tarball_and_loads_it(
    media_backend, tmp_path, mocker
):
    load = mocker.patch("viewer.management.commands.load_target.load_target")
    tarball = tmp_path / "Mpro.tgz"
    tarball.write_bytes(b"bundle")

    call_command("load_target", str(tarball), "-p", "lb12345-1")

    name = load.call_args.args[0]
    assert _upload_file_name(name) == "Mpro.tgz"
    assert media_storage.read_bytes(name) == b"bundle"
    assert load.call_args.kwargs["proposal_ref"] == "lb12345-1"


# Storing an upload's files -----------------------------------------------------


@pytest.fixture
def extracted_upload(tmp_path):
    """A loader (stand-in) with an extracted upload_2 and a stored bundle."""
    raw_data = tmp_path / "Mpro"
    (raw_data / "upload_2" / "aligned_files" / "x0001").mkdir(parents=True)
    (raw_data / "upload_2" / "aligned_files" / "x0001" / "x0001.pdb").write_text("A")
    (raw_data / "upload_2" / "meta_aligner.yaml").write_text("meta")

    def _make():
        stored_bundle = media_storage.store_bytes("tmp/uploads/u1/Mpro.tgz", b"bundle")
        return SimpleNamespace(
            raw_data=raw_data,
            version_dir="upload_2",
            final_path=Path("target_loader_data/Mpro_lb1"),
            data_bundle="Mpro.tgz",
            stored_bundle=stored_bundle,
            stored_prefix=None,
            final_bundle=None,
        )

    return _make


def test_store_upload_files_stores_the_upload_and_moves_the_bundle(
    media_backend, extracted_upload
):
    loader = extracted_upload()

    store_upload_files(loader)

    assert sorted(o.name for o in media_storage.list_prefix("target_loader_data")) == [
        "target_loader_data/Mpro_lb1/Mpro.tgz",
        "target_loader_data/Mpro_lb1/upload_2/aligned_files/x0001/x0001.pdb",
        "target_loader_data/Mpro_lb1/upload_2/meta_aligner.yaml",
    ]
    assert media_storage.read_bytes("target_loader_data/Mpro_lb1/Mpro.tgz") == b"bundle"
    assert not media_storage.exists(loader.stored_bundle)


def test_store_upload_files_does_not_replace_an_earlier_bundle(
    media_backend, extracted_upload
):
    media_storage.store_bytes("target_loader_data/Mpro_lb1/Mpro.tgz", b"first")
    loader = extracted_upload()

    store_upload_files(loader)

    assert media_storage.read_bytes("target_loader_data/Mpro_lb1/Mpro.tgz") == b"first"
    assert (
        media_storage.read_bytes("target_loader_data/Mpro_lb1/Mpro_upload_2.tgz")
        == b"bundle"
    )


def test_delete_stored_upload_files_deletes_only_this_upload(
    media_backend, extracted_upload
):
    media_storage.store_bytes("target_loader_data/Mpro_lb1/upload_1/a.pdb", b"1")
    loader = extracted_upload()
    store_upload_files(loader)

    delete_stored_upload_files(loader)

    assert [o.name for o in media_storage.list_prefix("target_loader_data")] == [
        "target_loader_data/Mpro_lb1/upload_1/a.pdb"
    ]


def test_delete_stored_upload_files_does_nothing_before_the_store(
    media_backend, extracted_upload
):
    media_storage.store_bytes("target_loader_data/Mpro_lb1/upload_2/a.pdb", b"x")
    loader = extracted_upload()

    delete_stored_upload_files(loader)

    assert media_storage.exists("target_loader_data/Mpro_lb1/upload_2/a.pdb")


# Explaining a failed store -------------------------------------------------------

_ABORT = (
    "ERROR: Aborting upload - failed to write"
    " 'target_loader_data/Mpro_lb1/upload_2/a+b.ccp4'"
    " to S3 bucket 'frag-media' (SignatureDoesNotMatch)"
)


def _s3_upload_error():
    return S3UploadError(
        "target_loader_data/Mpro_lb1/upload_2/a+b.ccp4",
        "frag-media",
        "SignatureDoesNotMatch",
    )


def test_a_failed_s3_store_is_explained_in_the_report():
    report = UploadReport(task=None, proposal_ref="lb1")

    report_store_failure(report, _s3_upload_error())

    assert report.json() == [_ABORT]
    assert report.failed


def test_other_failures_add_nothing_to_the_report():
    report = UploadReport(task=None, proposal_ref="lb1")

    report_store_failure(report, ValueError("bad data"))

    assert not report.json()


def _stored_bundle() -> str:
    """Stores a small bundle (one file in a Mpro/ directory). Returns its name."""
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode="w:gz") as tar:
        content = b"meta"
        info = tarfile.TarInfo("Mpro/upload_2/meta_aligner.yaml")
        info.size = len(content)
        tar.addfile(info, io.BytesIO(content))
    return media_storage.store_bytes("tmp/uploads/u1/Mpro.tgz", data.getvalue())


@pytest.mark.django_db
def test_load_target_reports_a_failed_s3_store(media_backend, mocker):
    mocker.patch.object(TargetLoader, "process_bundle")
    mocker.patch(
        "viewer.target_loader.store_upload_files", side_effect=_s3_upload_error()
    )
    final = mocker.patch.object(UploadReport, "final", autospec=True)

    load_target(_stored_bundle(), proposal_ref="lb1")

    report = final.call_args.args[0]
    assert _ABORT in report.json()
    assert final.call_args.kwargs == {"success": False}


# Reading .mol files ---------------------------------------------------------------


def test_read_stored_mol(media_backend):
    media_storage.store_bytes("computed_set_data/a.mol", _MOL_BLOCK)

    assert Chem.MolToSmiles(read_stored_mol("computed_set_data/a.mol")) == "CCO"
    assert read_stored_mol("computed_set_data/missing.mol") is None
    assert read_stored_mol(None) is None


@pytest.fixture
def loader(tmp_path):
    with mock.patch("viewer.target_loader.UploadReport"):
        target_loader = TargetLoader(str(tmp_path / "Mpro.tgz"), "lb12345-1", tmp_path)
    target_loader._final_path = Path(  # pylint: disable=protected-access
        "target_loader_data/Mpro_lb1"
    )
    return target_loader


def test_ligand_mol_is_read_from_the_extracted_upload(media_backend, loader):
    mol_file = loader.raw_data / "upload_2" / "x.mol"
    mol_file.parent.mkdir(parents=True)
    mol_file.write_text(_MOL_BLOCK)

    assert Chem.MolToSmiles(loader.read_ligand_mol("upload_2/x.mol")) == "CCO"


def test_ligand_mol_of_an_earlier_upload_is_read_from_storage(media_backend, loader):
    media_storage.store_bytes("target_loader_data/Mpro_lb1/upload_1/x.mol", _MOL_BLOCK)

    assert Chem.MolToSmiles(loader.read_ligand_mol("upload_1/x.mol")) == "CCO"


def test_a_missing_ligand_mol_is_none(media_backend, loader):
    assert loader.read_ligand_mol("upload_1/missing.mol") is None


# Target directory names -------------------------------------------------------------


def test_sanitize_directory_name_avoids_existing_names():
    assert sanitize_directory_name("Mpro lb1", existing=["Mpro_lb1"]) == "Mpro_lb1_2"
    assert sanitize_directory_name("Mpro lb1", existing=[]) == "Mpro_lb1"
