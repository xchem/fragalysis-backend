"""Tests for how computed-set uploads resolve ``ref_pdb`` (``MolOps``).

Without an uploaded PDB, ``ref_pdb`` names the observation whose protein the
virtual molecule sits in. It must resolve to exactly one current (not
superseded) observation of the target. A name that is ambiguous or matches
nothing fails the whole upload instead of guessing, or of silently creating a
molecule with no protein, crystal form site or conformer site.
"""
# pylint: disable=redefined-outer-name,unused-argument
from datetime import datetime, timezone

import pytest
from rdkit import Chem

from viewer.cset_upload import MolOps, blank_mol_vals
from viewer.models import (
    ComputedSet,
    Experiment,
    ExperimentUpload,
    SiteObservation,
    TagCategory,
)

# the categories TagManager.tag_new_site_observations looks up
TAG_CATEGORIES = (
    "Other",
    "CanonSites",
    "ConformerSites",
    "Quatassemblies",
    "Crystalforms",
    "CrystalformSites",
)

# ref_mols always resolves to this one, so only ref_pdb varies
INSPIRATION = "x0009a"

HEADER = {
    "ref_url": "https://example.com/method",
    "submitter_name": "Tester",
    "submitter_email": "tester@example.com",
    "submitter_institution": "Lab",
    "generation_date": "2026-10-29",
    "method": "refpdbtest",
}


@pytest.fixture
def ref_target(settings, tmp_path, user, make_project, make_target):
    """A target with one upload, and a factory adding observations to it."""
    media_root = tmp_path / "media"
    media_root.mkdir()
    settings.MEDIA_ROOT = str(media_root)

    project = make_project("proposal", members=[user])
    target = make_target(project, title="RefPdb")
    for category in TAG_CATEGORIES:
        TagCategory.objects.get_or_create(category=category)
    target.zip_archive = "RefPdb_proposal"
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

    def _observation(code, longcode=None, superseded=False):
        crystal = code.rstrip("abcdefghijklmnopqrstuvwxyz")
        experiment, _ = Experiment.objects.get_or_create(
            experiment_upload=upload,
            code=f"RefPdb-{crystal}",
            defaults={"pdb_info": f"crystallographic_files/{crystal}.pdb"},
        )
        return SiteObservation.objects.create(
            code=code,
            longcode=longcode,
            experiment=experiment,
            smiles="CCO",
            version=1,
            superseded=superseded,
        )

    _observation(INSPIRATION)
    target.observation = _observation
    return target


@pytest.fixture
def upload_ref_pdb(tmp_path, user, ref_target):
    """Upload an SDF with one data molecule per given ``ref_pdb`` value."""

    def _upload(*ref_pdbs):
        sdf = tmp_path / "upload.sdf"
        writer = Chem.SDWriter(str(sdf))
        header = Chem.MolFromSmiles("C")
        header.SetProp("_Name", "ver_1.2")
        for key, value in HEADER.items():
            header.SetProp(key, value)
        writer.write(header)

        for i, ref_pdb in enumerate(ref_pdbs):
            mol = Chem.MolFromSmiles("CC" + "C" * i)
            mol.SetProp("_Name", f"mol_{i}")
            mol.SetProp("ref_mols", INSPIRATION)
            mol.SetProp("ref_pdb", ref_pdb)
            writer.write(mol)
        writer.close()

        submitter_name, submitter_method, version = blank_mol_vals(str(sdf))
        computed_set, _ = MolOps(
            user_id=user.pk,
            sdf_filename=str(sdf),
            submitter_name=submitter_name,
            submitter_method=submitter_method,
            target=ref_target.pk,
            version=version,
            zfile=None,
            zfile_hashvals=None,
            computed_set_id=None,
        ).task()
        return computed_set

    return _upload


def _virtual_pdbs(computed_set):
    return [
        str(so.virtual_pdb_info)
        for so in SiteObservation.objects.filter(
            computedsetsiteobservation__computed_set=computed_set
        )
    ]


def _assert_nothing_uploaded():
    assert not ComputedSet.objects.exists()
    assert not SiteObservation.objects.filter(experiment__isnull=True).exists()


@pytest.mark.django_db
def test_unresolved_name_fails_the_whole_upload(ref_target, upload_ref_pdb):
    ref_target.observation("x0007a")

    with pytest.raises(ValueError, match="mol_1"):
        upload_ref_pdb("x0007a", "nosuch")

    _assert_nothing_uploaded()
