"""The media deletion watcher is quietened for deletions we asked for.

filewatcher.sh logs every file removed under media/target_loader_data, to catch
files vanishing when nothing should have touched them. A deliberate deletion -
delete_target, or peeling an upload off one - would bury that in thousands of
lines, so it drops a marker the watcher checks before writing.
"""

from pathlib import Path

from viewer.media_watcher import PAUSE_PREFIX, deletions_expected
from viewer.target_delete import delete_target

# pylint: disable=redefined-outer-name,unused-argument


def _markers(media_root) -> list:
    return sorted(Path(media_root).glob(f"{PAUSE_PREFIX}*"))


def test_a_marker_exists_only_for_the_duration(settings, tmp_path):
    settings.MEDIA_ROOT = str(tmp_path)

    assert not _markers(tmp_path)
    with deletions_expected("deleting something"):
        during = _markers(tmp_path)
        assert len(during) == 1
        assert "deleting something" in during[0].read_text()
    assert not _markers(tmp_path)


def test_the_marker_is_removed_even_when_the_deletion_fails(settings, tmp_path):
    """Otherwise a failed deletion would silence the watcher indefinitely."""
    settings.MEDIA_ROOT = str(tmp_path)

    try:
        with deletions_expected("deletion that explodes"):
            raise RuntimeError("boom")
    except RuntimeError:
        pass

    assert not _markers(tmp_path)


def test_the_marker_names_the_process(settings, tmp_path):
    """Concurrent deletions must not clear each other's marker."""
    import os

    settings.MEDIA_ROOT = str(tmp_path)

    with deletions_expected("mine"):
        assert _markers(tmp_path)[0].name.endswith(str(os.getpid()))


def test_an_unwritable_media_root_does_not_stop_the_deletion(settings, tmp_path):
    """Log tidiness must never be the reason a deletion fails."""
    settings.MEDIA_ROOT = str(tmp_path / "does" / "not" / "exist")

    with deletions_expected("nowhere to write a marker"):
        pass  # no exception is the assertion


def test_target_deletion_suppresses_the_watcher(
    db, settings, tmp_path, user, make_project, make_target, monkeypatch
):
    """delete_target runs inside the suppression, and leaves none behind."""
    media_root = tmp_path / "media"
    (media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY).mkdir(parents=True)
    (media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY).mkdir(parents=True)
    settings.MEDIA_ROOT = str(media_root)

    target = make_target(make_project("proposal", members=[user]), title="Quiet")
    target.zip_archive = "Quiet_proposal"
    target.save()
    subdir = media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY / "Quiet_proposal"
    subdir.mkdir()
    (subdir / "data.pdb").write_text("x")

    seen = []
    real_rmtree = __import__("shutil").rmtree

    def spy(path, *args, **kwargs):
        seen.append(_markers(media_root))
        return real_rmtree(path, *args, **kwargs)

    monkeypatch.setattr("viewer.target_delete.shutil.rmtree", spy)

    delete_target(target)

    # the marker was in place while the files were being removed ...
    assert seen and len(seen[0]) == 1
    # ... and is gone afterwards
    assert not _markers(media_root)
