"""Tests for the media-cleanup scan and the ``cleanup_media`` command.

The scan has to be exactly right in one direction: a file it fails to recognise
as referenced gets deleted. Most of what follows is therefore about the three
ways a computed-set reference can be stored (absolute path, ``computed_set_data/``
prefix, bare basename) and about the entries that must never be touched at all -
symlinks, the file watcher's log, anything a hand-edited ``zip_archive`` points at.

``django_cleanup`` is active (``fragalysis.settings``), so replacing a FileField
value also removes the old file from disk: tests that need a superseded file on
disk write it *after* the save.
"""

import os
from io import StringIO
from pathlib import Path

from django.core.management import call_command

from viewer.media_cleanup import (
    referenced_computed_set_names,
    scan_computed_set_data,
    scan_target_loader_data,
    unexpected_computed_set_fields,
)
from viewer.models import ComputedSet, SiteObservation, Target

# The `db` fixture is requested for its side effect (DB access) only, and
# fixtures legitimately reuse their names as arguments - both are standard
# pytest patterns that pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument


def _make_media_root(settings, tmp_path) -> Path:
    """Point MEDIA_ROOT at a temp dir with the standard media subdirectories."""
    media_root = tmp_path / "media"
    (media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY).mkdir(parents=True)
    (media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY).mkdir(parents=True)
    settings.MEDIA_ROOT = str(media_root)
    return media_root


def _loader_dir(settings, media_root) -> Path:
    return media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY


def _cset_dir(settings, media_root) -> Path:
    return media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY


def _target_with_dir(make_target, project, loader_dir, title, dir_name) -> Target:
    """A Target claiming ``dir_name``, with that directory laid down on disk."""
    target = make_target(project, title=title)
    target.zip_archive = dir_name
    target.save()
    subdir = loader_dir / dir_name
    subdir.mkdir()
    (subdir / "data.pdb").write_text("x")
    return target


def _orphan_dir(loader_dir, name="Gone_proposal") -> Path:
    """A target-loader directory no Target claims."""
    path = loader_dir / name
    path.mkdir()
    (path / "data.pdb").write_text("y")
    return path


def _run(*argv) -> str:
    """Run the command, returning everything it wrote.

    Always pass string arguments: ``call_command`` with keywords bypasses
    argparse, and the defaults with it.
    """
    out = StringIO()
    call_command("cleanup_media", *argv, stdout=out, stderr=out)
    return out.getvalue()


# --------------------------------------------------------------------------- #
# Reporting and deleting
# --------------------------------------------------------------------------- #


def test_report_only_lists_orphans_and_deletes_nothing(
    db, settings, tmp_path, user, make_project, make_target
):
    """The default run is a report: it names the orphans but writes nothing."""
    media_root = _make_media_root(settings, tmp_path)
    project = make_project("proposal", members=[user])
    _target_with_dir(
        make_target, project, _loader_dir(settings, media_root), "Keep", "Keep_proposal"
    )
    orphan_dir = _orphan_dir(_loader_dir(settings, media_root))
    orphan_file = _cset_dir(settings, media_root) / "orphan.sdf"
    orphan_file.write_text("orphan")

    output = _run()

    assert orphan_dir.exists()
    assert orphan_file.exists()
    assert "Gone_proposal" in output
    assert "orphan.sdf" in output
    assert "--delete" in output


def test_delete_removes_orphans_and_keeps_referenced(
    db, settings, tmp_path, user, make_project, make_target
):
    """--delete removes the unreferenced directory and leaves the claimed one."""
    media_root = _make_media_root(settings, tmp_path)
    loader_dir = _loader_dir(settings, media_root)
    project = make_project("proposal", members=[user])
    target = _target_with_dir(make_target, project, loader_dir, "Keep", "Keep_proposal")
    orphan_dir = _orphan_dir(loader_dir)

    _run("--delete")

    assert not orphan_dir.exists()
    assert (loader_dir / "Keep_proposal" / "data.pdb").exists()
    assert Target.objects.filter(pk=target.pk).exists()


def test_deletions_log_is_never_removed(
    db, settings, tmp_path, user, make_project, make_target
):
    """filewatcher.sh writes its log into target_loader_data/ - a file there is
    reported, never deleted, however orphaned it looks."""
    media_root = _make_media_root(settings, tmp_path)
    loader_dir = _loader_dir(settings, media_root)
    project = make_project("proposal", members=[user])
    _target_with_dir(make_target, project, loader_dir, "Keep", "Keep_proposal")
    watcher_log = loader_dir / "deletions-blackadder.log"
    watcher_log.write_text("[2026-08-07] Media watcher started")

    output = _run("--delete")

    assert watcher_log.exists()
    assert "deletions-blackadder.log" in output
    assert "not a directory" in output


# --------------------------------------------------------------------------- #
# Computed-set reference resolution - the part that must not get this wrong
# --------------------------------------------------------------------------- #


def test_every_reference_flavour_keeps_its_file(db, settings, tmp_path, user):
    """All three storage flavours protect their file, and only their file.

    The migrated ``virtual_ligand_mol`` rows (bare basename, file in
    computed_set_data/) are the ones a naive "resolve under MEDIA_ROOT" rule
    would wrongly reclaim.
    """
    media_root = _make_media_root(settings, tmp_path)
    cset_dir = _cset_dir(settings, media_root)

    # (a) written_sdf_filename, absolute - plus the bare submitted_sdf naming
    # the *original* upload, whose real file carries the mangled name.
    written = cset_dir / "setA_upload_1_orig.sdf"
    written.write_text("A")
    ComputedSet.objects.create(
        name="setA",
        md_ordinal=1,
        submitted_sdf="orig.sdf",
        written_sdf_filename=str(written),
        owner_user=user,
    )

    # (b) the same, with no written_sdf_filename: only the reconstruction from
    # (name, md_ordinal, submitted_sdf) can protect this one.
    reconstructed = cset_dir / "setB_upload_2_orig2.sdf"
    reconstructed.write_text("B")
    ComputedSet.objects.create(
        name="setB", md_ordinal=2, submitted_sdf="orig2.sdf", owner_user=user
    )

    # (c) an assay-data set: a properly saved computed_set_data/ value, and no
    # name or ordinal to reconstruct from.
    assay = cset_dir / "assay.csv"
    assay.write_text("C")
    ComputedSet.objects.create(
        submitted_sdf=f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/assay.csv",
        owner_user=user,
    )

    # (d) virtual_pdb_info, computed_set_data/-prefixed.
    pdb = cset_dir / "A0486a#deadbeef.pdb_cafe"
    pdb.write_text("D")
    SiteObservation.objects.create(
        virtual_pdb_info=f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/{pdb.name}"
    )

    # (e) virtual_ligand_mol as migration 0149 wrote it: a bare basename.
    migrated = cset_dir / "setA_upload_1_x_y_z.mol"
    migrated.write_text("E")
    SiteObservation.objects.create(virtual_ligand_mol=migrated.name)

    # (f) virtual_ligand_mol as cset_upload writes it today: the file lives in
    # the target-loader tree, so it protects nothing here - a same-named file in
    # computed_set_data/ is still debris.
    same_name_elsewhere = cset_dir / "current.mol"
    same_name_elsewhere.write_text("F")
    SiteObservation.objects.create(
        virtual_ligand_mol=(
            f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/T/upload_1/virtual_files/"
            "current.mol"
        )
    )

    orphan = cset_dir / "nobody_points_here.sdf"
    orphan.write_text("G")

    _run("--delete")

    for survivor in (written, reconstructed, assay, pdb, migrated):
        assert survivor.exists(), f"{survivor.name} was wrongly deleted"
    assert not same_name_elsewhere.exists()
    assert not orphan.exists()


def test_missing_referenced_file_is_reported_not_fatal(db, settings, tmp_path, user):
    """A reference with no file on disk shows up as 'referenced but absent'."""
    media_root = _make_media_root(settings, tmp_path)
    ComputedSet.objects.create(
        name="setA",
        md_ordinal=1,
        submitted_sdf=(f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/never_written.sdf"),
        owner_user=user,
    )
    orphan = _cset_dir(settings, media_root) / "orphan.sdf"
    orphan.write_text("x")

    report = scan_computed_set_data()

    assert "never_written.sdf" in report.missing
    assert [c.path for c in report.candidates] == [orphan]


def test_superseded_value_becomes_debris(db, settings, tmp_path, user):
    """Replacing a field's value leaves the old file unreferenced, so the next
    run reclaims it while the current one survives."""
    media_root = _make_media_root(settings, tmp_path)
    cset_dir = _cset_dir(settings, media_root)
    prefix = settings.COMPUTED_SET_MEDIA_DIRECTORY

    computed_set = ComputedSet.objects.create(
        submitted_sdf=f"{prefix}/first.sdf", owner_user=user
    )
    computed_set.submitted_sdf = f"{prefix}/second.sdf"
    computed_set.save()

    # Written after the save: django_cleanup removes the replaced file.
    superseded = cset_dir / "first.sdf"
    superseded.write_text("old")
    (cset_dir / "second.sdf").write_text("new")

    assert "first.sdf" not in referenced_computed_set_names()

    _run("--delete")

    assert not superseded.exists()
    assert (cset_dir / "second.sdf").exists()


def test_no_unexpected_computed_set_fields():
    """Tripwire: if this fails, a new FileField points at computed_set_data/ and
    its files would be treated as debris. Add it to ``_COMPUTED_SET_FIELDS`` in
    viewer/media_cleanup.py."""
    assert unexpected_computed_set_fields() == []


# --------------------------------------------------------------------------- #
# Safety rules
# --------------------------------------------------------------------------- #


def test_empty_database_refuses_to_delete_without_force(db, settings, tmp_path):
    """No targets at all is the signature of the wrong database, not of 25 GB of
    debris - so it takes --force."""
    media_root = _make_media_root(settings, tmp_path)
    orphan_dir = _orphan_dir(_loader_dir(settings, media_root))
    orphan_file = _cset_dir(settings, media_root) / "orphan.sdf"
    orphan_file.write_text("x")

    output = _run("--delete")

    assert orphan_dir.exists()
    assert orphan_file.exists()
    assert "Refusing to delete" in output

    _run("--delete", "--force")

    assert not orphan_dir.exists()
    assert not orphan_file.exists()


def test_symlinks_are_never_followed_or_deleted(db, settings, tmp_path):
    """Download staging trees are built from symlinks; a symlink is reported and
    left, and whatever it points at is untouched."""
    media_root = _make_media_root(settings, tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    precious = outside / "precious.pdb"
    precious.write_text("do not delete")

    dir_link = _loader_dir(settings, media_root) / "Linked_proposal"
    dir_link.symlink_to(outside, target_is_directory=True)
    file_link = _cset_dir(settings, media_root) / "linked.sdf"
    file_link.symlink_to(precious)

    output = _run("--delete", "--force")

    assert precious.exists()
    assert dir_link.is_symlink()
    assert file_link.is_symlink()
    assert "symlink" in output


def test_unusable_zip_archive_is_ignored_with_a_warning(
    db, settings, tmp_path, user, make_project, make_target
):
    """A zip_archive that isn't a plain directory name never matches anything."""
    media_root = _make_media_root(settings, tmp_path)
    loader_dir = _loader_dir(settings, media_root)
    project = make_project("proposal", members=[user])

    traversal = make_target(project, title="Traversal")
    traversal.zip_archive = "../outside"
    traversal.save()
    legacy = make_target(project, title="Legacy")
    legacy.zip_archive = "None"
    legacy.save()

    sibling = media_root / "outside"
    sibling.mkdir()
    (sibling / "keep.txt").write_text("keep")
    orphan_dir = _orphan_dir(loader_dir)

    report = scan_target_loader_data()
    _run("--delete")

    assert sibling.exists() and (sibling / "keep.txt").exists()
    assert not orphan_dir.exists()
    assert any("unusable zip_archive" in w for w in report.warnings)
    assert any("empty zip_archive" in w for w in report.warnings)


def test_rows_referencing_files_inside_hold_the_directory(
    db, settings, tmp_path, user, make_project, make_target
):
    """A directory no Target claims but whose files rows still use is an
    inconsistency to report, not something to delete."""
    media_root = _make_media_root(settings, tmp_path)
    loader_dir = _loader_dir(settings, media_root)
    project = make_project("proposal", members=[user])
    make_target(project, title="Mpro")

    stranded = _orphan_dir(loader_dir, name="Stranded_proposal")
    SiteObservation.objects.create(
        bound_file=(
            f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/Stranded_proposal/data.pdb"
        )
    )

    output = _run("--delete")

    assert stranded.exists()
    assert "rows reference files inside it" in output


def test_scope_flags_limit_the_area(db, settings, tmp_path):
    """Each area can be cleaned on its own."""
    media_root = _make_media_root(settings, tmp_path)
    orphan_dir = _orphan_dir(_loader_dir(settings, media_root))
    orphan_file = _cset_dir(settings, media_root) / "orphan.sdf"
    orphan_file.write_text("x")

    _run("--delete", "--force", "--computed-set-data")
    assert orphan_dir.exists()
    assert not orphan_file.exists()

    _run("--delete", "--force", "--target-loader-data")
    assert not orphan_dir.exists()


def test_missing_media_directories_are_a_noop(db, settings, tmp_path):
    """A MEDIA_ROOT without either subdirectory scans clean rather than raising."""
    settings.MEDIA_ROOT = str(tmp_path / "empty-media")

    output = _run("--delete", "--force")

    assert scan_target_loader_data().candidates == []
    assert scan_computed_set_data().candidates == []
    assert "removed 0 entries" in output


def test_directory_size_ignores_symlinked_content(db, settings, tmp_path):
    """A symlink inside an orphan contributes its own size, not its target's."""
    media_root = _make_media_root(settings, tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    big = outside / "big.bin"
    big.write_bytes(b"0" * 100_000)

    orphan = _orphan_dir(_loader_dir(settings, media_root))
    os.symlink(big, orphan / "link.bin")

    candidates = scan_target_loader_data().candidates

    assert len(candidates) == 1
    assert candidates[0].size_bytes < 1_000
    assert big.exists()
