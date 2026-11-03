# pylint: disable=redefined-outer-name,unused-argument,protected-access
"""Tests for target and upload deletes through the media storage seam (issue 1038).

Each S3 request takes 30-75 ms, so files are deleted in batches:

- django_cleanup does not delete SiteObservation and ExperimentUpload files one at
  a time (``@cleanup.ignore``). Observations are deleted only with their target
  (or, rarely, with an upload), and those deletes remove the files themselves.
- A target delete removes ``target_loader_data/<target>/`` with ``delete_prefix``
  and the target's computed-set files with one ``delete_many``.
- An upload delete removes its directory's unreferenced files with ``delete_many``,
  and its own bundle (never an earlier upload's bundle with the same name).

Each test runs with both media storages (the ``media_backend`` fixture).
"""
from datetime import datetime, timezone

import pytest
from django.core.files.storage import storages

from viewer import media_storage
from viewer.models import (
    ComputedSet,
    ComputedSetSiteObservation,
    Experiment,
    ExperimentUpload,
    SiteObservation,
)
from viewer.target_delete import delete_target
from viewer.upload_delete import (
    NotDeletable,
    UploadDeletionPlan,
    _bundle_archive_name,
    _remove_upload_media,
    read_upload_metadata,
)

_DIR = "target_loader_data/Mpro_lb1"


@pytest.fixture
def target(user, make_project, make_target):
    target = make_target(make_project("lb1", members=[user]), title="Mpro")
    target.zip_archive = "Mpro_lb1"
    target.save()
    return target


def _upload(target, user, upload_data_dir="upload_1", version=1, file="Mpro.tgz"):
    return ExperimentUpload.objects.create(
        project=target.project,
        target=target,
        file=file,
        upload_data_dir=upload_data_dir,
        upload_version=version,
        commit_datetime=datetime(2026, 1, version, tzinfo=timezone.utc),
        committer=user,
    )


def _all_names() -> list[str]:
    return sorted(o.name for o in media_storage.list_prefix(""))


# django_cleanup ------------------------------------------------------------------------


def test_deleting_an_observation_keeps_its_files(media_backend, target, user):
    name = media_storage.store_bytes(f"{_DIR}/upload_1/x0001_apo.pdb", b"A")
    obs = SiteObservation.objects.create(apo_file=name)

    obs.delete()

    assert media_storage.exists(name)


def test_deleting_an_upload_keeps_its_files(media_backend, target, user):
    name = media_storage.store_bytes(f"{_DIR}/upload_1/assemblies.yaml", b"A")
    upload = _upload(target, user)
    upload.assembly_transforms = name
    upload.save()

    upload.delete()

    assert media_storage.exists(name)


# Target delete -----------------------------------------------------------------------------


@pytest.fixture
def loaded_target(target, user, settings):
    """A target with observations and a computed set, all with stored files."""
    upload = _upload(target, user)
    experiment = Experiment.objects.create(experiment_upload=upload, code="x0001")
    for n in range(3):
        files = {
            field: media_storage.store_bytes(f"{_DIR}/upload_1/x{n}_{field}.pdb", b"x")
            for field in ("apo_file", "bound_file", "apo_desolv_file")
        }
        SiteObservation.objects.create(experiment=experiment, code=f"x{n}", **files)

    cset_dir = settings.COMPUTED_SET_MEDIA_DIRECTORY
    written = media_storage.store_bytes(f"{cset_dir}/cset_upload_1_hits.sdf", b"s")
    computed_set = ComputedSet.objects.create(
        name="cset",
        md_ordinal=1,
        target=target,
        submitted_sdf="hits.sdf",
        written_sdf_filename=written,
        owner_user=user,
    )
    pdb = media_storage.store_bytes(f"{cset_dir}/cset_mol.pdb", b"p")
    obs = SiteObservation.objects.create(virtual_pdb_info=pdb)
    ComputedSetSiteObservation.objects.create(
        computed_set=computed_set, site_observation=obs
    )

    # Another target's files (in a sibling directory and in computed_set_data/)
    media_storage.store_bytes("target_loader_data/Mpro_lb1_2/upload_1/y.pdb", b"y")
    media_storage.store_bytes(f"{cset_dir}/other_upload_1_hits.sdf", b"o")
    return target


def test_delete_target_removes_its_files_only(media_backend, loaded_target, settings):
    delete_target(loaded_target)

    cset_dir = settings.COMPUTED_SET_MEDIA_DIRECTORY
    assert _all_names() == [
        f"{cset_dir}/other_upload_1_hits.sdf",
        "target_loader_data/Mpro_lb1_2/upload_1/y.pdb",
    ]


@pytest.mark.parametrize("media_backend", ["s3"], indirect=True)
def test_delete_target_deletes_in_batches(media_backend, loaded_target, mocker):
    client = storages["default"].client
    single = mocker.spy(client, "delete_object")
    batches = mocker.spy(client, "delete_objects")

    delete_target(loaded_target)

    # 12 files: not one request each. (The Target's own few FileFields may
    # still be deleted singly by django_cleanup.)
    assert single.call_count <= 2
    assert batches.call_count >= 1


# Upload delete ----------------------------------------------------------------------------------


def test_read_upload_metadata_reads_the_stored_file(
    media_backend, target, user, settings
):
    upload = _upload(target, user)
    media_storage.store_bytes(
        f"{_DIR}/upload_1/meta_aligner.yaml",
        f"data_format_version: {settings.XCA_DATA_FORMAT_VERSION}\n"
        "version_number: 1\ncrystals: {}\n",
    )

    assert isinstance(read_upload_metadata(upload), dict)


def test_read_upload_metadata_without_the_file_is_not_deletable(
    media_backend, target, user
):
    with pytest.raises(NotDeletable):
        read_upload_metadata(_upload(target, user))


def test_remove_upload_media_keeps_referenced_files(
    media_backend, target, user, mocker
):
    upload = _upload(target, user, upload_data_dir="upload_2", version=2)
    kept = media_storage.store_bytes(f"{_DIR}/upload_2/kept.pdb", b"k")
    gone = media_storage.store_bytes(f"{_DIR}/upload_2/sub/gone.pdb", b"g")
    media_storage.store_bytes(f"{_DIR}/upload_1/earlier.pdb", b"e")
    mocker.patch("viewer.upload_delete._referenced_names", return_value={kept})
    plan = UploadDeletionPlan(upload.pk, upload.upload_version, "Mpro")

    _remove_upload_media(upload, plan)

    assert _all_names() == [f"{_DIR}/upload_1/earlier.pdb", kept]
    assert plan.files_removed == [gone]
    assert plan.files_kept == [kept]
    assert any("could not be removed" in w for w in plan.warnings)


def test_remove_upload_media_removes_the_whole_directory(
    media_backend, target, user, mocker
):
    upload = _upload(target, user, upload_data_dir="upload_2", version=2)
    media_storage.store_bytes(f"{_DIR}/upload_2/a.pdb", b"a")
    media_storage.store_bytes(f"{_DIR}/upload_2/sub/b.pdb", b"b")
    mocker.patch("viewer.upload_delete._referenced_names", return_value=set())
    plan = UploadDeletionPlan(upload.pk, upload.upload_version, "Mpro")

    _remove_upload_media(upload, plan)

    assert not list(media_storage.list_prefix(f"{_DIR}/upload_2"))
    assert media_storage.list_dirs(_DIR) == []
    assert plan.warnings == []


def test_the_bundle_of_an_upload_is_its_own(media_backend, target, user):
    media_storage.store_bytes(f"{_DIR}/Mpro.tgz", b"1")
    media_storage.store_bytes(f"{_DIR}/Mpro_upload_2.tgz", b"2")
    first = _upload(target, user)
    second = _upload(target, user, upload_data_dir="upload_2", version=2)

    assert _bundle_archive_name(second) == f"{_DIR}/Mpro_upload_2.tgz"
    assert _bundle_archive_name(first) == f"{_DIR}/Mpro.tgz"


def test_an_earlier_uploads_bundle_is_not_deleted(media_backend, target, user):
    """upload_2's renamed bundle is missing, so "Mpro.tgz" is upload_1's."""
    media_storage.store_bytes(f"{_DIR}/Mpro.tgz", b"1")
    _upload(target, user)
    second = _upload(target, user, upload_data_dir="upload_2", version=2)

    assert _bundle_archive_name(second) is None


def test_remove_upload_media_removes_the_bundle(media_backend, target, user, mocker):
    media_storage.store_bytes(f"{_DIR}/Mpro.tgz", b"1")
    media_storage.store_bytes(f"{_DIR}/Mpro_upload_2.tgz", b"2")
    _upload(target, user)
    second = _upload(target, user, upload_data_dir="upload_2", version=2)
    mocker.patch("viewer.upload_delete._referenced_names", return_value=set())
    plan = UploadDeletionPlan(second.pk, second.upload_version, "Mpro")

    _remove_upload_media(second, plan)

    assert _all_names() == [f"{_DIR}/Mpro.tgz"]
    assert plan.files_removed == [f"{_DIR}/Mpro_upload_2.tgz"]
