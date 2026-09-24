"""An earlier upload's crystallographic files must survive a later upload.

django_cleanup deletes the file behind a FileField the moment the field's value is
replaced. A re-upload that re-supplies a crystal rewrites Experiment.pdb_info (and
mtz_info/cif_info) from the earlier upload's directory to its own, which destroyed the
earlier upload's copy on disk while its meta_aligner.yaml still listed it.

Found by loading the six real A71EV2A bundles: A71EV2A-x0152.pdb is in the upload_1
archive and on disk after uploads 1 and 2, and gone after upload 3 - the first upload to
change that crystal's pdb path. The event maps beside it survived only because map_info
is an ArrayField, which django_cleanup does not track.

Guarded by @cleanup.ignore on Experiment; these tests fail without it.
"""

from pathlib import Path

import pytest

from viewer.models import Experiment, ExperimentUpload

# django_cleanup defers its deletion to transaction.on_commit, so the ordinary `db`
# fixture - which rolls back - never runs it and these tests would pass either way.
# `django_capture_on_commit_callbacks(execute=True)` runs those callbacks without a real
# commit. Do NOT reach for `transactional_db` here: it truncates every table afterwards,
# including rows seeded by data migrations (service_status.ServiceState), and with
# --reuse-db they do not come back - which breaks unrelated tests for the rest of the run.
# pylint: disable=redefined-outer-name,unused-argument


@pytest.fixture
def experiment_with_file(db, settings, tmp_path, user, make_project, make_target):
    """An Experiment whose pdb_info points at a real file, as a load leaves it."""
    settings.MEDIA_ROOT = str(tmp_path)
    project = make_project("proposal", members=[user])
    target = make_target(project, title="Keeper")
    upload = ExperimentUpload.objects.create(
        project=project,
        target=target,
        committer=user,
        commit_datetime="2026-01-01T00:00:00Z",
        upload_data_dir="upload_1",
        upload_version=1,
    )
    relative = "target_loader_data/Keeper/upload_1/crystallographic_files/x0152.pdb"
    absolute = Path(tmp_path) / relative
    absolute.parent.mkdir(parents=True)
    absolute.write_text("the original crystal file")

    experiment = Experiment.objects.create(
        experiment_upload=upload, code="x0152", pdb_info=relative
    )
    return experiment, absolute


def test_replacing_a_file_path_keeps_the_previous_file(
    experiment_with_file, django_capture_on_commit_callbacks
):
    """The case that bit us: a later upload re-supplies the crystal."""
    experiment, original = experiment_with_file

    with django_capture_on_commit_callbacks(execute=True):
        experiment.pdb_info = (
            "target_loader_data/Keeper/upload_3/crystallographic_files/x0152.pdb"
        )
        experiment.save(update_fields=["pdb_info"])

    assert (
        original.is_file()
    ), "upload_1's copy was deleted when upload_3 re-supplied it"
    assert original.read_text() == "the original crystal file"


def test_deleting_the_experiment_keeps_the_file(
    experiment_with_file, django_capture_on_commit_callbacks
):
    """Removing the file is the caller's job now - target_delete and upload_delete.

    Pinned so nobody reintroduces automatic deletion by a different route: those two
    know which upload a file belongs to, and this signal does not.
    """
    experiment, original = experiment_with_file

    with django_capture_on_commit_callbacks(execute=True):
        experiment.delete()

    assert original.is_file()
