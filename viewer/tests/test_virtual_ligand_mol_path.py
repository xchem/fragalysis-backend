"""Tests for resolving ``SiteObservation.virtual_ligand_mol`` to a real path.

The field holds three historical flavours of value (see
``viewer/media_cleanup.py`` for the same catalogue):

- a **bare basename**, written by migration 0149 when ComputedMolecule was
  merged into SiteObservation, with the file itself in ``computed_set_data/``;
- a **MEDIA_ROOT-relative path**, which is what ``viewer.cset_upload`` writes
  today - ``target_loader_data/<archive>/upload_N/virtual_files/<name>.mol``
  (assigning the name directly bypasses the field's ``upload_to``);
- an **absolute path**, for completeness.

``SiteObservation.get_virtual_ligand_mol_path()`` is the single resolver both
``viewer.download_structures`` and ``viewer.target_loader`` use. Before it
existed, the download joined every value onto ``computed_set_data/`` and so
raised ``OSError: Bad input file media/computed_set_data/target_loader_data/...``
for anything cset_upload had written - see issue #1025.
"""

from pathlib import Path
from types import SimpleNamespace

from rdkit import Chem

from viewer.download_structures import DownloadStructures
from viewer.models import ComputedSet, ComputedSetSiteObservation, SiteObservation

# The `db` fixture is requested for its side effect (DB access) only, and
# fixtures legitimately reuse their names as arguments - both standard pytest
# patterns that pylint misreads. _compound_sets_zip is the unit under test, so
# reaching for it directly is deliberate.
# pylint: disable=redefined-outer-name,unused-argument,protected-access

MOL_BLOCK = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))


def test_resolves_bare_basename_into_computed_set_data(settings, tmp_path):
    """Migration 0149 rows: the file lives in ``computed_set_data/``."""
    settings.MEDIA_ROOT = str(tmp_path)
    obs = SiteObservation(virtual_ligand_mol="A0486a.mol")

    assert obs.get_virtual_ligand_mol_path() == Path(
        tmp_path, settings.COMPUTED_SET_MEDIA_DIRECTORY, "A0486a.mol"
    )


def test_resolves_media_relative_path_against_media_root(settings, tmp_path):
    """What cset_upload writes today: a path relative to MEDIA_ROOT, in the
    target-loader tree - it must *not* be joined onto ``computed_set_data/``."""
    settings.MEDIA_ROOT = str(tmp_path)
    name = (
        f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/A71EV2A/upload_1/virtual_files/x.mol"
    )
    obs = SiteObservation(virtual_ligand_mol=name)

    assert obs.get_virtual_ligand_mol_path() == Path(tmp_path, name)


def test_resolves_computed_set_relative_path_against_media_root(settings, tmp_path):
    """A ``computed_set_data/...`` value is also MEDIA_ROOT-relative, so it
    resolves once - not twice (which would double the directory)."""
    settings.MEDIA_ROOT = str(tmp_path)
    name = f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/A0486a.mol"
    obs = SiteObservation(virtual_ligand_mol=name)

    assert obs.get_virtual_ligand_mol_path() == Path(tmp_path, name)


def test_resolves_absolute_path_unchanged(settings, tmp_path):
    settings.MEDIA_ROOT = str(tmp_path)
    obs = SiteObservation(virtual_ligand_mol="/elsewhere/A0486a.mol")

    assert obs.get_virtual_ligand_mol_path() == Path("/elsewhere/A0486a.mol")


def test_returns_none_when_unset(settings, tmp_path):
    """An unset FileField stringifies to '' (and, historically, to 'None')."""
    settings.MEDIA_ROOT = str(tmp_path)

    assert SiteObservation().get_virtual_ligand_mol_path() is None
    assert SiteObservation(virtual_ligand_mol="").get_virtual_ligand_mol_path() is None
    assert (
        SiteObservation(virtual_ligand_mol="None").get_virtual_ligand_mol_path() is None
    )


def _downloader(target, tempdir):
    """A DownloadStructures with just enough of a Celery task to log with."""
    return DownloadStructures(
        tempdir=str(tempdir),
        task=SimpleNamespace(request=SimpleNamespace(id="test-task")),
        target=target,
        use_zip=False,
        target_access_string="proposal",
    )


def _computed_set(user, target, observations):
    cset = ComputedSet.objects.create(
        name="cset", target=target, owner_user=user, submitted_sdf="hits.sdf"
    )
    for obs in observations:
        ComputedSetSiteObservation.objects.create(
            computed_set=cset, site_observation=obs
        )
    return cset


def test_compound_sets_zip_reads_cset_upload_written_mol(
    settings, tmp_path, user, make_project, make_target
):
    """The regression from #1025: a mol written by cset_upload (and so stored
    under ``target_loader_data/``) must be found and written to the archive."""
    media_root = tmp_path / "media"
    virtual_dir = (
        media_root
        / settings.TARGET_LOADER_MEDIA_DIRECTORY
        / "A71EV2A"
        / "upload_1"
        / "virtual_files"
    )
    virtual_dir.mkdir(parents=True)
    (virtual_dir / "x.mol").write_text(MOL_BLOCK)
    settings.MEDIA_ROOT = str(media_root)

    target = make_target(make_project("proposal", members=[user]), title="A71EV2A")
    obs = SiteObservation.objects.create(
        smiles="CCO",
        virtual_name="hit-1",
        virtual_ligand_mol=str(
            Path(
                settings.TARGET_LOADER_MEDIA_DIRECTORY,
                "A71EV2A",
                "upload_1",
                "virtual_files",
                "x.mol",
            )
        ),
    )
    _computed_set(user, target, [obs])

    tempdir = tmp_path / "download"
    tempdir.mkdir()
    _downloader(target, tempdir)._compound_sets_zip(target)

    written = (tempdir / "virtual_hits" / "hits.sdf").read_text()
    assert "hit-1" in written


def test_compound_sets_zip_skips_unreadable_mol(
    settings, tmp_path, user, make_project, make_target
):
    """A missing or unparseable mol is skipped, not fatal: the remaining
    observations still make it into the archive."""
    media_root = tmp_path / "media"
    computed_dir = media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY
    computed_dir.mkdir(parents=True)
    (computed_dir / "good.mol").write_text(MOL_BLOCK)
    (computed_dir / "corrupt.mol").write_text("not a molfile")
    settings.MEDIA_ROOT = str(media_root)

    target = make_target(make_project("proposal", members=[user]), title="A71EV2A")
    # Bare basenames, as migration 0149 wrote them.
    good = SiteObservation.objects.create(
        smiles="CCO", virtual_name="good-hit", virtual_ligand_mol="good.mol"
    )
    corrupt = SiteObservation.objects.create(
        smiles="CCO", virtual_name="corrupt-hit", virtual_ligand_mol="corrupt.mol"
    )
    missing = SiteObservation.objects.create(
        smiles="CCO", virtual_name="missing-hit", virtual_ligand_mol="gone.mol"
    )
    unset = SiteObservation.objects.create(smiles="CCO", virtual_name="unset-hit")
    _computed_set(user, target, [good, corrupt, missing, unset])

    tempdir = tmp_path / "download"
    tempdir.mkdir()
    _downloader(target, tempdir)._compound_sets_zip(target)

    written = (tempdir / "virtual_hits" / "hits.sdf").read_text()
    assert "good-hit" in written
    for skipped in ("corrupt-hit", "missing-hit", "unset-hit"):
        assert skipped not in written
