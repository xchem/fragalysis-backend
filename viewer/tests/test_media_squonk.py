# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the Squonk2 file transfer and job-file upload (issue 1038).

- A file transfer gets the requested media files into scratch (in parallel),
  then gives DmApi the local copies.
- A job-file upload pulls the job's SDF and .meta.json from Squonk into
  scratch, adds the blank molecule and the annotations, then stores the SDF as
  ``tmp/squonk_upload/<code>/job_job-<code>.sdf``. The Celery chain passes that
  storage name to validate_compound_set, and the upload's files are deleted
  at the end.

DmApi and the database records are mocked. Each test runs with both media
storages (the ``media_backend`` fixture).
"""
import json
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest

from viewer import media_storage
from viewer.squonk_job_file_transfer import process_file_transfer
from viewer.squonk_job_file_upload import (
    delete_upload_files,
    get_upload_sub_directory,
    process_compound_set_file,
)
from viewer.utils import SDF_VERSION

_CODE = "1b2c3d4e-0000-0000-0000-000000000001"
_SDF = "mol-1\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n"
_JOB_END = datetime(2026, 11, 2, 10, 0, tzinfo=timezone.utc)
_META = {"annotations": [{"fields": {"DDG": {"description": "Delta G"}}}]}


# File transfer -----------------------------------------------------------------------


@pytest.fixture
def transfer(mocker):
    """A mocked JobFileTransfer and DmApi. Returns (record, put_files mock, copies)."""
    record = SimpleNamespace(
        target="Mpro",
        squonk_project="project-1",
        sub_path="abcd",
        proteins=[],
        compounds=[],
        transfer_progress=0,
        save=lambda: None,
    )
    mocker.patch(
        "viewer.squonk_job_file_transfer.JobFileTransfer.objects.get",
        return_value=record,
    )
    copies: dict[str, bytes] = {}

    def _put(auth_token, project_id, project_files, project_path, force):
        # Read the local copies while they exist (they are scratch files)
        for path in project_files:
            copies[Path(path).name] = Path(path).read_bytes()
        return SimpleNamespace(success=True, msg="")

    put = mocker.patch(
        "viewer.squonk_job_file_transfer.DmApi.put_unmanaged_project_files",
        side_effect=_put,
    )
    return record, put, copies


def test_transfer_gives_dmapi_local_copies(media_backend, transfer):
    record, put, copies = transfer
    media_storage.store_bytes("target_loader_data/Mpro/a/x0001_apo-desolv.pdb", b"P")
    media_storage.store_bytes("computed_set_data/v1a b#1.mol", b"M")
    record.proteins = ["target_loader_data/Mpro/a/x0001_apo-desolv.pdb"]
    record.compounds = ["computed_set_data/v1a%20b%231.mol"]  # URL-encoded

    process_file_transfer("token", 1)

    assert copies == {"x0001_apo-desolv.pdb": b"P", "v1a b#1.mol": b"M"}
    assert put.call_args.kwargs["project_path"] == "/fragalysis-files/abcd"
    assert record.transfer_progress == 100
    # The scratch copies are removed
    assert not any(Path(p).exists() for p in put.call_args.kwargs["project_files"])


def test_transfer_of_a_missing_file_fails(media_backend, transfer):
    record, put, _ = transfer
    record.proteins = ["target_loader_data/Mpro/a/missing.pdb"]

    with pytest.raises(RuntimeError, match="No such file"):
        process_file_transfer("token", 1)

    assert put.call_count == 0


# Job-file upload ---------------------------------------------------------------------------


@pytest.fixture
def job_request(mocker):
    record = SimpleNamespace(
        id=1,
        code=_CODE,
        squonk_project="project-1",
        squonk_url_ext="data-manager-ui/results/instance/i-1",
        squonk_job_info={"msg": {"callback_token": "callback"}},
        user=SimpleNamespace(username="alan", email="alan@example.com"),
    )
    mocker.patch(
        "viewer.squonk_job_file_upload.JobRequest.objects.get", return_value=record
    )
    mocker.patch(
        "viewer.squonk_job_file_upload._SQ2A.get_ui_url",
        return_value="https://squonk.example.com/",
    )
    return record


@pytest.fixture
def squonk_files(mocker):
    """DmApi gives the job's files. Set an entry to None to make it fail."""
    files = {"out.meta.json": json.dumps(_META), "out.sdf": _SDF}

    def _get(token, project_id, project_path, project_file, local_file, timeout_s):
        content = files.get(project_file)
        if content is None:
            return SimpleNamespace(success=False)
        Path(local_file).write_text(content, encoding="utf-8")
        return SimpleNamespace(success=True)

    mocker.patch(
        "viewer.squonk_job_file_upload.DmApi.get_unmanaged_project_file_with_token",
        side_effect=_get,
    )
    return files


def test_upload_stores_the_processed_sdf(media_backend, job_request, squonk_files):
    name = process_compound_set_file(1, _JOB_END, "/out", "out.sdf")

    assert name == f"tmp/squonk_upload/{_CODE}/job_job-{_CODE}.sdf"
    sdf = media_storage.read_text(name)
    assert sdf.startswith(SDF_VERSION)
    assert "<DDG>" in sdf
    assert "mol-1" in sdf
    # Only the processed SDF is stored
    assert [o.name for o in media_storage.list_prefix("tmp")] == [name]


@pytest.mark.parametrize("missing", ["out.meta.json", "out.sdf"])
def test_upload_with_a_missing_job_file_stores_nothing(
    media_backend, job_request, squonk_files, missing
):
    squonk_files[missing] = None

    name = process_compound_set_file(1, _JOB_END, "/out", "out.sdf")

    assert name == f"tmp/squonk_upload/{_CODE}/job_job-{_CODE}.sdf"
    assert not media_storage.exists(name)


def test_delete_upload_files_deletes_only_this_job(media_backend, job_request):
    media_storage.store_bytes(f"{get_upload_sub_directory(job_request)}/a.sdf", b"x")
    media_storage.store_bytes("tmp/squonk_upload/other/a.sdf", b"keep")

    delete_upload_files(job_request)

    assert [o.name for o in media_storage.list_prefix("tmp")] == [
        "tmp/squonk_upload/other/a.sdf"
    ]
