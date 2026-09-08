"""Tests pinning the one path convention for ``virtual_ligand_mol``.

``SiteObservation.virtual_ligand_mol`` holds a path relative to MEDIA_ROOT.
One convention, resolved one way, by everybody.

It did not always. Migration 0149 wrote a bare basename with the file in
``computed_set_data/``; commit 726624a4 then moved the files into the
target-loader tree and switched the column to a MEDIA_ROOT-relative path
without updating either reader or backfilling the migrated rows. The column
ended up holding two incompatible shapes, and
``viewer.download_structures._compound_sets_zip`` guessed wrong for everything
``viewer.cset_upload`` had written - ``OSError: Bad input file
media/computed_set_data/target_loader_data/...``, issue #1025.

The fix normalises the data (migration 0172) rather than teaching callers to
detect which shape they are holding. These tests exist to keep it that way:
``virtual_ligand_mol_path`` must stay a plain MEDIA_ROOT join, because the
moment it starts interpreting the value the two conventions become permanent.
"""

import importlib
from pathlib import Path
from types import SimpleNamespace

import pytest
from rdkit import Chem

from viewer.download_structures import DownloadStructures
from viewer.models import ComputedSet, ComputedSetSiteObservation, SiteObservation

# The `db` fixture is requested for its side effect (DB access) only, and
# fixtures legitimately reuse their names as arguments - both standard pytest
# patterns that pylint misreads. _compound_sets_zip is the unit under test, so
# reaching for it directly is deliberate.
# pylint: disable=redefined-outer-name,unused-argument,protected-access

MOL_BLOCK = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))


# --------------------------------------------------------------------------
# the accessor
# --------------------------------------------------------------------------


def test_resolves_media_relative_value(settings, tmp_path):
    """The convention: MEDIA_ROOT + the stored value, and nothing else."""
    settings.MEDIA_ROOT = str(tmp_path)
    name = (
        f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/A71EV2A/upload_1/virtual_files/x.mol"
    )

    obs = SiteObservation(virtual_ligand_mol=name)

    assert obs.virtual_ligand_mol_path == Path(tmp_path, name)


def test_computed_set_value_resolves_once_not_twice(settings, tmp_path):
    """A ``computed_set_data/...`` value is MEDIA_ROOT-relative like any other -
    it must not be joined onto ``computed_set_data/`` a second time."""
    settings.MEDIA_ROOT = str(tmp_path)
    name = f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/A0486a.mol"

    obs = SiteObservation(virtual_ligand_mol=name)

    assert obs.virtual_ligand_mol_path == Path(tmp_path, name)


def test_returns_none_when_unset(settings, tmp_path):
    settings.MEDIA_ROOT = str(tmp_path)

    assert SiteObservation().virtual_ligand_mol_path is None
    assert SiteObservation(virtual_ligand_mol="").virtual_ligand_mol_path is None


def test_accessor_does_not_rescue_a_bare_basename(settings, tmp_path):
    """The guard against reintroducing flavour-detection.

    A bare basename is not a valid value any more - migration 0172 removed
    them. If this ever starts resolving into ``computed_set_data/``, the
    accessor has gone back to guessing and the column is free to hold two
    conventions again.
    """
    settings.MEDIA_ROOT = str(tmp_path)

    obs = SiteObservation(virtual_ligand_mol="A0486a.mol")

    assert obs.virtual_ligand_mol_path == Path(tmp_path, "A0486a.mol")


# --------------------------------------------------------------------------
# migration 0172
# --------------------------------------------------------------------------


@pytest.fixture
def run_normalisation():
    """Call migration 0172's data function against the real models.

    The function only uses ``get_model``, so a stub ``apps`` is enough and we
    avoid dragging the whole migration graph into a unit test.
    """
    # import_module, not `from ... import`: the module name starts with a
    # digit, so it has no importable dotted form.
    migration = importlib.import_module(
        "viewer.migrations.0172_normalise_virtual_ligand_mol"
    )

    def _run():
        migration.normalise_virtual_ligand_mol(
            SimpleNamespace(get_model=lambda *_: SiteObservation), None
        )

    return _run


def test_migration_normalises_bare_basenames(db, run_normalisation, settings):
    """The 3817-row case: a basename becomes ``computed_set_data/<name>``,
    which points at the file migration 0149 actually wrote."""
    obs = SiteObservation.objects.create(smiles="CCO", virtual_ligand_mol="A0486a.mol")

    run_normalisation()

    obs.refresh_from_db()
    assert (
        obs.virtual_ligand_mol.name
        == f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/A0486a.mol"
    )


def test_migration_leaves_current_values_alone(db, settings, run_normalisation):
    """Values already in the convention are untouched - the Tibor-stack case."""
    name = (
        f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/A71EV2A/upload_1/virtual_files/x.mol"
    )
    obs = SiteObservation.objects.create(smiles="CCO", virtual_ligand_mol=name)
    unset = SiteObservation.objects.create(smiles="CCC")

    run_normalisation()

    obs.refresh_from_db()
    unset.refresh_from_db()
    assert obs.virtual_ligand_mol.name == name
    assert not unset.virtual_ligand_mol


def test_migration_is_idempotent(db, settings, run_normalisation):
    """Running it twice must not produce ``computed_set_data/computed_set_data/``."""
    SiteObservation.objects.create(smiles="CCO", virtual_ligand_mol="A0486a.mol")

    run_normalisation()
    run_normalisation()

    obs = SiteObservation.objects.get(smiles="CCO")
    assert (
        obs.virtual_ligand_mol.name
        == f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/A0486a.mol"
    )


# --------------------------------------------------------------------------
# the download itself
# --------------------------------------------------------------------------


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


def test_compound_sets_zip_reads_a_cset_upload_written_mol(
    settings, tmp_path, user, make_project, make_target
):
    """The #1025 regression: a mol written by cset_upload, and so stored under
    ``target_loader_data/``, must be found and written to the archive."""
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


def test_compound_sets_zip_skips_unusable_mols(
    settings, tmp_path, user, make_project, make_target
):
    """A missing or unparseable mol is skipped, not fatal: the remaining
    observations still make it into the archive.

    Before this, a missing file raised OSError out of MolFromMolFile and an
    unparseable one gave `mol = None` and then AttributeError on SetProp -
    either way the whole download died.
    """
    media_root = tmp_path / "media"
    computed_dir = media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY
    computed_dir.mkdir(parents=True)
    (computed_dir / "good.mol").write_text(MOL_BLOCK)
    (computed_dir / "corrupt.mol").write_text("not a molfile")
    settings.MEDIA_ROOT = str(media_root)

    target = make_target(make_project("proposal", members=[user]), title="A71EV2A")
    prefix = settings.COMPUTED_SET_MEDIA_DIRECTORY
    good = SiteObservation.objects.create(
        smiles="CCO", virtual_name="good-hit", virtual_ligand_mol=f"{prefix}/good.mol"
    )
    corrupt = SiteObservation.objects.create(
        smiles="CCO",
        virtual_name="corrupt-hit",
        virtual_ligand_mol=f"{prefix}/corrupt.mol",
    )
    missing = SiteObservation.objects.create(
        smiles="CCO",
        virtual_name="missing-hit",
        virtual_ligand_mol=f"{prefix}/gone.mol",
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
