"""Tests for the media files a computed-set upload writes (``MolOps``).

The per-molecule ``.mol`` files and the PDBs from an uploaded zip are staged
and only moved into MEDIA_ROOT once the upload's transaction has committed. A
failed upload, whatever the exception, leaves neither them nor the extracted
zip PDBs behind.
"""
# pylint: disable=redefined-outer-name,unused-argument
from datetime import datetime, timezone

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from viewer import cset_upload
from viewer.cset_upload import MolOps, blank_mol_vals
from viewer.models import Experiment, ExperimentUpload, SiteObservation, TagCategory

# the categories TagManager.tag_new_site_observations looks up
TAG_CATEGORIES = (
    "Other",
    "CanonSites",
    "ConformerSites",
    "Quatassemblies",
    "Crystalforms",
    "CrystalformSites",
)

REF_CODE = "x0001a"

HEADER = {
    "ref_url": "https://example.com/method",
    "submitter_name": "Tester",
    "submitter_email": "tester@example.com",
    "submitter_institution": "Lab",
    "generation_date": "2026-10-29",
    "method": "mediatest",
}


@pytest.fixture
def media_root(settings, tmp_path):
    root = tmp_path / "media"
    root.mkdir()
    settings.MEDIA_ROOT = str(root)
    return root


@pytest.fixture
def media_target(media_root, user, make_project, make_target):
    """A target with one upload and one observation to reference."""
    project = make_project("proposal", members=[user])
    target = make_target(project, title="Media")
    for category in TAG_CATEGORIES:
        TagCategory.objects.get_or_create(category=category)
    target.zip_archive = "Media_proposal"
    target.save()

    upload = ExperimentUpload.objects.create(
        project=project,
        target=target,
        committer=user,
        commit_datetime=datetime(2026, 1, 1, tzinfo=timezone.utc),
        upload_data_dir="upload_1",
        upload_version=1,
        data_version_major=3,
        data_version_minor=1,
        file="bundle_upload_1.tgz",
    )
    experiment = Experiment.objects.create(
        experiment_upload=upload,
        code="Media-x0001",
        pdb_info="crystallographic_files/x0001.pdb",
    )
    SiteObservation.objects.create(
        code=REF_CODE, experiment=experiment, smiles="CCO", version=1
    )
    return target


@pytest.fixture
def upload_mols(tmp_path, user, media_target):
    """Upload an SDF with one data molecule per (smiles, ref_pdb) pair."""

    def _upload(*mols, zfile=None):
        sdf = tmp_path / "upload.sdf"
        writer = Chem.SDWriter(str(sdf))
        header = Chem.MolFromSmiles("C")
        header.SetProp("_Name", "ver_1.2")
        for key, value in HEADER.items():
            header.SetProp(key, value)
        writer.write(header)

        for i, (smiles, ref_pdb) in enumerate(mols):
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            AllChem.EmbedMolecule(mol, randomSeed=1)
            mol = Chem.RemoveHs(mol)
            mol.SetProp("_Name", f"mol_{i}")
            mol.SetProp("ref_mols", REF_CODE)
            mol.SetProp("ref_pdb", ref_pdb)
            writer.write(mol)
        writer.close()

        submitter_name, submitter_method, version = blank_mol_vals(str(sdf))
        computed_set, _ = MolOps(
            user_id=user.pk,
            sdf_filename=str(sdf),
            submitter_name=submitter_name,
            submitter_method=submitter_method,
            target=media_target.pk,
            version=version,
            zfile=zfile,
            zfile_hashvals={} if zfile is not None else None,
            computed_set_id=None,
        ).task()
        return computed_set

    return _upload


@pytest.fixture
def zip_pdb(media_root):
    """A PDB as PdbOps.save_pdb_zip leaves it: in MEDIA_ROOT/tmp, keyed by code."""
    (media_root / "tmp").mkdir(exist_ok=True)
    (media_root / "tmp" / "upload1.pdb").write_text("ATOM\n", encoding="utf-8")
    return {"upload1": "tmp/upload1.pdb"}


def _files(media_root, subdir):
    return sorted(p for p in (media_root / subdir).rglob("*") if p.is_file())


def _assert_no_staging_left(media_root):
    tmp = media_root / "tmp"
    assert not tmp.exists() or not any(p.is_dir() for p in tmp.iterdir())


@pytest.mark.django_db
def test_mol_files_are_moved_into_place(media_root, upload_mols):
    computed_set = upload_mols(("CCN", REF_CODE), ("CCC", REF_CODE))

    observations = SiteObservation.objects.filter(
        computedsetsiteobservation__computed_set=computed_set
    )
    written = {media_root / str(so.virtual_ligand_mol) for so in observations}
    assert len(written) == 2
    assert all(p.is_file() for p in written)
    assert set(_files(media_root, "target_loader_data")) == written
    _assert_no_staging_left(media_root)


@pytest.mark.django_db
def test_second_pose_of_a_compound_reads_the_first_staged_file(media_root, upload_mols):
    computed_set = upload_mols(("CCN", REF_CODE), ("CCN", REF_CODE))

    observations = SiteObservation.objects.filter(
        computedsetsiteobservation__computed_set=computed_set
    )
    assert observations.count() == 2
    assert len(_files(media_root, "target_loader_data")) == 2


@pytest.mark.django_db
def test_failed_upload_leaves_no_mol_files(media_root, upload_mols):
    with pytest.raises(ValueError, match="nosuch"):
        upload_mols(("CCN", REF_CODE), ("CCC", "nosuch"))

    assert not _files(media_root, "target_loader_data")
    _assert_no_staging_left(media_root)


@pytest.mark.django_db
def test_unexpected_exception_leaves_no_mol_files(media_root, upload_mols, monkeypatch):
    def _boom(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(cset_upload.TagManager, "tag_new_site_observations", _boom)

    with pytest.raises(RuntimeError, match="boom"):
        upload_mols(("CCN", REF_CODE))

    assert not _files(media_root, "target_loader_data")
    _assert_no_staging_left(media_root)


@pytest.mark.django_db
def test_zip_pdb_is_moved_into_place(media_root, upload_mols, zip_pdb):
    computed_set = upload_mols(("CCN", "upload1"), zfile=zip_pdb)

    so = SiteObservation.objects.get(
        computedsetsiteobservation__computed_set=computed_set
    )
    assert str(so.virtual_pdb_info).startswith("computed_set_data/upload1.pdb_")
    assert (media_root / str(so.virtual_pdb_info)).is_file()
    assert not (media_root / "tmp" / "upload1.pdb").exists()


@pytest.mark.django_db
def test_failed_upload_removes_zip_pdbs_and_keeps_its_message(
    media_root, upload_mols, zip_pdb
):
    with pytest.raises(ValueError, match="nosuch"):
        upload_mols(("CCN", "upload1"), ("CCC", "nosuch"), zfile=zip_pdb)

    assert not (media_root / "tmp" / "upload1.pdb").exists()
    assert not (media_root / "computed_set_data").exists() or not _files(
        media_root, "computed_set_data"
    )
    _assert_no_staging_left(media_root)
