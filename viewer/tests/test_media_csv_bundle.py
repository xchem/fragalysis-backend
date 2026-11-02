# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the DictToCSV and experiment-bundle downloads (issue 1038).

- DictToCSV stores its CSV as ``downloads/dicttocsv/<uuid>/<file>``. The GET
  returns the file only if the name is such a file, then deletes it.
- The experiment-bundle download finds the bundle where the loader stores it:
  ``target_loader_data/<target>/<file>`` (or ``<stem>_<upload_N><suffix>``,
  when a bundle of that name was already there), and nginx serves it.

Each test runs with both media storages (the ``media_backend`` fixture).
"""
import re
from datetime import datetime, timezone

import pytest

from viewer import media_storage
from viewer.models import ExperimentUpload

_CSV_NAME = re.compile(r"^downloads/dicttocsv/[0-9a-f-]{36}/(?P<file>[^/]+)$")


def _csv_file_name(name: str) -> str:
    match = _CSV_NAME.match(name)
    assert match, name
    return match["file"]


# DictToCSV ------------------------------------------------------------------------------


def _post_csv(client, filename="sites.csv"):
    response = client.post(
        "/api/dicttocsv/",
        {"title": "Sites", "filename": filename, "dict": {"a": [1, 2], "b": [3, 4]}},
        format="json",
    )
    assert response.status_code == 200
    return response.json()["file_url"]


def test_the_csv_is_stored_under_dicttocsv(media_backend, authenticated_client):
    file_url = _post_csv(authenticated_client)

    assert _csv_file_name(file_url) == "sites.csv"
    assert media_storage.read_text(file_url) == "Sites\na,b\n1,3\n2,4\n"


@pytest.mark.parametrize("filename", ["../../evil.csv", "/etc/evil.csv"])
def test_the_csv_file_name_stays_in_its_directory(
    media_backend, authenticated_client, filename
):
    file_url = _post_csv(authenticated_client, filename=filename)

    assert _csv_file_name(file_url) == "evil.csv"


@pytest.mark.parametrize("legacy_absolute", [False, True])
def test_get_returns_the_csv_then_deletes_it(
    media_backend, settings, authenticated_client, legacy_absolute
):
    file_url = _post_csv(authenticated_client)
    requested = (
        f"{settings.MEDIA_ROOT.rstrip('/')}/{file_url}" if legacy_absolute else file_url
    )

    response = authenticated_client.get(f"/api/dicttocsv/?file_url={requested}")

    assert response.status_code == 200
    assert response.content == b"Sites\na,b\n1,3\n2,4\n"
    assert response["Content-Disposition"] == "attachment; filename=sites.csv"
    assert list(media_storage.list_prefix("downloads/dicttocsv")) == []


@pytest.mark.parametrize(
    "file_url",
    [
        "downloads/t1/Mpro.tar.gz",
        "downloads/dicttocsv/x/../../t1/Mpro.tar.gz",
        "downloads/dicttocsv",
        "../secret",
    ],
)
def test_get_refuses_other_files(media_backend, authenticated_client, file_url):
    media_storage.store_bytes("downloads/t1/Mpro.tar.gz", b"archive")

    response = authenticated_client.get(f"/api/dicttocsv/?file_url={file_url}")

    assert response.json() == "Please provide a file_url for an existing DictToCsv file"
    assert media_storage.exists("downloads/t1/Mpro.tar.gz")


def test_get_of_a_missing_csv(media_backend, authenticated_client):
    response = authenticated_client.get(
        "/api/dicttocsv/?file_url=downloads/dicttocsv/"
        "12345678-1234-1234-1234-123456789abc/sites.csv"
    )

    assert response.json() == "The given DictToCsv file does not exist"


# The experiment bundle ------------------------------------------------------------------


@pytest.fixture
def target(user, make_project, make_target):
    target = make_target(make_project("lb1", members=[user]), title="Mpro")
    target.zip_archive = "Mpro_lb1"
    target.save()
    return target


def _upload(target, user, upload_data_dir="upload_1", file="Mpro.tgz"):
    return ExperimentUpload.objects.create(
        project=target.project,
        target=target,
        file=file,
        upload_data_dir=upload_data_dir,
        commit_datetime=datetime(2026, 1, 1, tzinfo=timezone.utc),
        committer=user,
    )


def test_bundle_name_is_in_the_target_directory(media_backend, target, user):
    upload = _upload(target, user)
    media_storage.store_bytes("target_loader_data/Mpro_lb1/Mpro.tgz", b"1")

    assert upload.bundle_name() == "target_loader_data/Mpro_lb1/Mpro.tgz"


def test_bundle_name_of_a_second_upload_with_the_same_name(media_backend, target, user):
    first = _upload(target, user)
    second = _upload(target, user, upload_data_dir="upload_2")
    media_storage.store_bytes("target_loader_data/Mpro_lb1/Mpro.tgz", b"1")
    media_storage.store_bytes("target_loader_data/Mpro_lb1/Mpro_upload_2.tgz", b"2")

    assert first.bundle_name() == "target_loader_data/Mpro_lb1/Mpro.tgz"
    assert second.bundle_name() == "target_loader_data/Mpro_lb1/Mpro_upload_2.tgz"


def test_bundle_name_of_a_missing_bundle_is_none(media_backend, target, user):
    assert _upload(target, user).bundle_name() is None


def _post_bundle(client, target, filename="Mpro.tgz"):
    return client.post(
        "/api/download_target_experiments/",
        {"project": target.project.pk, "target": target.pk, "filename": filename},
        format="json",
    )


def test_the_bundle_is_served_through_nginx(
    media_backend, target, user, authenticated_client
):
    _upload(target, user)
    media_storage.store_bytes("target_loader_data/Mpro_lb1/Mpro.tgz", b"bundle")

    response = _post_bundle(authenticated_client, target)

    assert response.status_code == 200
    location = "/media-s3/" if media_backend == "s3" else "/media-local/"
    assert response["X-Accel-Redirect"] == (
        location + "target_loader_data/Mpro_lb1/Mpro.tgz"
    )
    assert response["Content-Disposition"] == 'attachment; filename="Mpro.tgz"'


def test_a_missing_bundle_is_not_found(
    media_backend, target, user, authenticated_client
):
    _upload(target, user)

    response = _post_bundle(authenticated_client, target)

    assert response.status_code == 404


def test_a_target_in_another_project_is_refused(
    media_backend, target, user, make_project, authenticated_client
):
    other = make_project("lb2", members=[user])
    _upload(target, user)

    response = authenticated_client.post(
        "/api/download_target_experiments/",
        {"project": other.pk, "target": target.pk, "filename": "Mpro.tgz"},
        format="json",
    )

    assert response.status_code == 403
    assert response.json() == {"error": "The Target is not part of the given Project"}
