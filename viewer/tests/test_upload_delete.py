"""Tests for the upload-deletion feature.

Covers the reusable service (``viewer.upload_delete.delete_latest_upload``) and the
``delete_target`` management command whose default is to peel off the latest upload.

The fixtures are hand-built rather than loaded through the target loader:
``test_target_loader.py`` is skipped wholesale for want of a committed bundle, so nothing
here can be driven through a real load. Each graph is written together with the
``meta_aligner.yaml`` that would have produced it, because that file - not the database -
is what the deleter reads to decide what the upload added.

The test settings leave ``DEPLOYMENT_MODE`` at DEVELOPMENT, so the command's production
guard is off unless a test turns it on.
"""

from datetime import datetime, timezone
from pathlib import Path

import pytest
import yaml
from django.core.management import call_command
from django.core.management.base import CommandError

from viewer.models import (
    CanonSite,
    CanonSiteConf,
    Compound,
    Experiment,
    ExperimentCompound,
    ExperimentUpload,
    Pose,
    SiteObservation,
    Xtalform,
    XtalformSite,
)
from viewer.upload_delete import (
    NotDeletable,
    delete_latest_upload,
    plan_upload_deletion,
    read_upload_metadata,
    restore_superseded,
)

# The `db` fixture is requested for its side effect (DB access) only, and fixtures
# legitimately reuse their names as arguments - both are standard pytest patterns that
# pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument

ALTLOC = "0"
CHAIN = "A"
LIGAND = "501"


# --------------------------------------------------------------------------- #
# Fixture construction
# --------------------------------------------------------------------------- #


def _longcode(code: str, version: int) -> str:
    return f"{code}_{CHAIN}_{LIGAND}_{ALTLOC}_v{version}"


def _write_metadata(
    directory: Path,
    version_number: int,
    observations: list[tuple[str, int]],
    canon_sites: list[tuple[str, int]],
    conformer_sites: list[tuple[str, int]],
    xtalform_sites: list[tuple[str, int]],
    data_format_version: str = "3.1",
    crystal_files: dict | None = None,
) -> Path:
    """Write a trimmed dfv 3.1 ``meta_aligner.yaml`` describing ``observations``.

    ``observations`` are ``(experiment code, version)`` pairs, nested the way a real
    bundle nests them: ``crystal -> aligned_files -> chain -> ligand -> altloc ->
    version``. The metadata is cumulative, so it lists carried-over entries at their
    original version alongside the ones this upload created.
    """
    crystals: dict = {}
    for code, version in observations:
        versions = (
            crystals.setdefault(code, {"aligned_files": {}})["aligned_files"]
            .setdefault(CHAIN, {})
            .setdefault(LIGAND, {})
            .setdefault(ALTLOC, {})
        )
        versions[version] = {"structure": f"{code}/{version}.pdb"}

    # Which upload currently supplies each crystal's own files. A re-upload that
    # re-supplies a carried-over crystal rewrites the Experiment's path to its own
    # directory, and this is the record of what it was before that.
    for code, supplied_by in (crystal_files or {}).items():
        crystals.setdefault(code, {"aligned_files": {}})["crystallographic_files"] = {
            "xtal_pdb": {
                "file": f"upload_{supplied_by}/crystallographic_files/{code}/{code}.pdb"
            }
        }

    def section(entries, separator):
        return {f"{name}{separator}{version}": {} for name, version in entries}

    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "meta_aligner.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "data_format_version": data_format_version,
                "version_number": version_number,
                "crystals": crystals,
                "canon_sites": section(canon_sites, "+"),
                "conformer_sites": section(conformer_sites, "+"),
                "xtalform_sites": section(xtalform_sites, "/"),
            }
        ),
        encoding="utf-8",
    )
    return path


@pytest.fixture
def make_graph(db, settings, tmp_path, user, make_project, make_target):
    """Factory building a multi-upload target, its media tree and its metadata.

    Returns a simple namespace-ish dict so tests can reach the individual rows. The graph
    is the smallest one that exercises every hazard the deleter has to handle:

    * an experiment carried over from upload 1 that gains a *new version* of its
      observation in upload 2 (so the old row must survive and be un-superseded);
    * an experiment created only by upload 2 (so it must go);
    * a pose whose main is the upload-2 row and whose members include the upload-1 row;
    * a pose whose only member is an upload-2 row (so it empties out);
    * a surviving CanonSiteConf whose ``ref_site_observation`` is a doomed row;
    * a CanonSiteConf created by upload 2 that a *surviving* observation was re-linked to;
    * media under each upload directory, including a file in upload 2 that a surviving
      experiment still points at.
    """

    def _make(uploads=2):
        media_root = tmp_path / "media"
        loader_dir = media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY
        loader_dir.mkdir(parents=True)
        settings.MEDIA_ROOT = str(media_root)

        project = make_project("proposal", members=[user])
        target = make_target(project, title="Peel")
        target.zip_archive = "Peel_proposal"
        target.save()
        target_dir = loader_dir / "Peel_proposal"

        xtalform = Xtalform.objects.create(name="xtalform1")
        canon_site = CanonSite.objects.create(
            name="cs1", residues=[], version=1, canon_site_num=1
        )
        xtalform_site = XtalformSite.objects.create(
            xtalform_site_id="xs1",
            xtalform=xtalform,
            canon_site=canon_site,
            lig_chain=CHAIN,
            residues=[],
            version=1,
        )
        conf_site = CanonSiteConf.objects.create(
            name="csc1", canon_site=canon_site, residues=[], version=1
        )
        compound = Compound.objects.create(
            smiles="CCO", inchi="InChI=1S/C2H6O", project=project
        )

        experiment_uploads = []
        for n in range(1, uploads + 1):
            upload_dir = target_dir / f"upload_{n}"
            upload_dir.mkdir(parents=True)
            experiment_uploads.append(
                ExperimentUpload.objects.create(
                    project=project,
                    target=target,
                    committer=user,
                    commit_datetime=datetime(2026, 1, n, tzinfo=timezone.utc),
                    upload_data_dir=f"upload_{n}",
                    upload_version=n,
                    data_version_major=3,
                    data_version_minor=1,
                    # The loader stores only the basename here, while it moves the
                    # archive into the target's directory - see _bundle_archive_paths.
                    file=f"bundle_upload_{n}.tgz",
                )
            )

        # The carried-over experiment, created by upload 1.
        carried = Experiment.objects.create(
            experiment_upload=experiment_uploads[0],
            code="Peel-x0001",
            xtalform=xtalform,
        )
        pose = Pose.objects.create(
            canon_site=canon_site, compound=compound, display_name="Peel-x0001_A_501"
        )

        observations = {}
        for n in range(1, uploads + 1):
            rel = f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/Peel_proposal/upload_{n}"
            bound = f"{rel}/Peel-x0001_v{n}.pdb"
            (media_root / bound).write_text(f"bound {n}")
            observations[n] = SiteObservation.objects.create(
                code=f"Peel-x0001_v{n}",
                longcode=_longcode("Peel-x0001", n),
                experiment=carried,
                cmpd=compound,
                xtalform_site=xtalform_site,
                canon_site_conf=conf_site,
                pose=pose,
                seq_id=int(LIGAND),
                chain_id=CHAIN,
                altloc=ALTLOC,
                smiles="CCO",
                version=n,
                # every version but the last is superseded
                superseded=n < uploads,
                bound_file=bound,
            )
        pose.main_site_observation = observations[uploads]
        pose.save()

        # An experiment that only the latest upload knows about, with its own pose.
        latest = experiment_uploads[-1]
        newcomer = Experiment.objects.create(
            experiment_upload=latest,
            code="Peel-x0002",
            xtalform=xtalform,
        )
        rel = (
            f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/Peel_proposal/"
            f"upload_{uploads}"
        )
        bound = f"{rel}/Peel-x0002_v{uploads}.pdb"
        (media_root / bound).write_text("bound newcomer")
        doomed_pose = Pose.objects.create(
            canon_site=canon_site, compound=compound, display_name="Peel-x0002_A_501"
        )
        newcomer_obs = SiteObservation.objects.create(
            code=f"Peel-x0002_v{uploads}",
            longcode=_longcode("Peel-x0002", uploads),
            experiment=newcomer,
            cmpd=compound,
            xtalform_site=xtalform_site,
            canon_site_conf=conf_site,
            pose=doomed_pose,
            seq_id=int(LIGAND),
            chain_id=CHAIN,
            altloc=ALTLOC,
            smiles="CCO",
            version=uploads,
            bound_file=bound,
        )
        doomed_pose.main_site_observation = newcomer_obs
        doomed_pose.save()

        # A *new version* of the conformer site, created by the latest upload, which the
        # surviving upload-1 observation was re-linked to - the loader recomputes this FK
        # on every load. A successor shares its predecessor's name and bumps the version,
        # the way the bundle keys them ("<name>+<version>").
        new_conf_site = CanonSiteConf.objects.create(
            name="csc1", canon_site=canon_site, residues=[], version=uploads
        )
        observations[1].canon_site_conf = new_conf_site
        observations[1].save(update_fields=["canon_site_conf"])

        # A surviving parent whose reference FK points at a doomed row.
        conf_site.ref_site_observation = observations[uploads]
        conf_site.save(update_fields=["ref_site_observation"])

        # The carried-over crystal's own files, re-supplied by every upload - which is
        # what rewrites the Experiment's path into the newest upload's directory.
        for n in range(1, uploads + 1):
            crystal_dir = (
                target_dir / f"upload_{n}" / "crystallographic_files" / "Peel-x0001"
            )
            crystal_dir.mkdir(parents=True, exist_ok=True)
            (crystal_dir / "Peel-x0001.pdb").write_text(f"crystal file from upload {n}")
        carried.pdb_info = (
            f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/Peel_proposal/"
            f"upload_{uploads}/crystallographic_files/Peel-x0001/Peel-x0001.pdb"
        )
        carried.save(update_fields=["pdb_info"])

        # Debris in the latest upload that nothing references: this must go.
        (media_root / f"{rel}/unreferenced.map").write_text("debris")

        # The bundle archives, where _move_and_save_target_experiment puts them: in the
        # target directory, NOT where ExperimentUpload.file's basename resolves to.
        for n in range(1, uploads + 1):
            (target_dir / f"bundle_upload_{n}.tgz").write_bytes(b"archive")

        # Metadata, written as the loader would have left it. Cumulative: every upload's
        # file lists everything the target holds at that point.
        for n in range(1, uploads + 1):
            obs_entries = [("Peel-x0001", v) for v in range(1, n + 1)]
            conf_entries = [("csc1", 1)]
            if n == uploads and uploads > 1:
                obs_entries.append(("Peel-x0002", uploads))
                conf_entries.append(("csc1", uploads))
            # Each upload re-supplies the carried-over crystal, so its metadata names
            # its own directory; the newcomer arrives with the last upload.
            supplied = {"Peel-x0001": n}
            if n == uploads and uploads > 1:
                supplied["Peel-x0002"] = uploads
            _write_metadata(
                target_dir / f"upload_{n}",
                version_number=n,
                observations=obs_entries,
                canon_sites=[("cs1", 1)],
                conformer_sites=conf_entries,
                xtalform_sites=[("xs1", 1)],
                crystal_files=supplied,
            )

        return {
            "media_root": media_root,
            "target_dir": target_dir,
            "target": target,
            "project": project,
            "uploads": experiment_uploads,
            "carried": carried,
            "newcomer": newcomer,
            "observations": observations,
            "newcomer_obs": newcomer_obs,
            "pose": pose,
            "doomed_pose": doomed_pose,
            "canon_site": canon_site,
            "conf_site": conf_site,
            "new_conf_site": new_conf_site,
            "xtalform_site": xtalform_site,
            "compound": compound,
        }

    return _make


# --------------------------------------------------------------------------- #
# Reading the bundle's metadata
# --------------------------------------------------------------------------- #


def test_read_upload_metadata_returns_only_what_the_upload_created(make_graph):
    """The file is cumulative; only entries at this upload's version were created here."""
    graph = make_graph()
    declared = read_upload_metadata(graph["uploads"][1])

    assert declared["SiteObservation"] == {
        _longcode("Peel-x0001", 2),
        _longcode("Peel-x0002", 2),
    }
    # Declared at version 1, so carried over rather than created by upload 2.
    assert _longcode("Peel-x0001", 1) not in declared["SiteObservation"]
    assert declared["CanonSiteConf"] == {"csc1"}
    assert declared["CanonSite"] == set()
    assert declared["XtalformSite"] == set()


def test_missing_metadata_is_refused(make_graph):
    """The bundle's metadata is the authority, so its absence is fatal, not a fallback."""
    graph = make_graph()
    (graph["target_dir"] / "upload_2" / "meta_aligner.yaml").unlink()

    with pytest.raises(NotDeletable, match="cannot establish what this upload added"):
        delete_latest_upload(graph["uploads"][1])


def test_old_data_format_is_refused(make_graph):
    """dfv 2.2 nests aligned_files without an altloc level; refuse rather than misparse."""
    graph = make_graph()
    _write_metadata(
        graph["target_dir"] / "upload_2",
        version_number=2,
        observations=[("Peel-x0001", 2)],
        canon_sites=[],
        conformer_sites=[],
        xtalform_sites=[],
        data_format_version="2.2",
    )

    with pytest.raises(NotDeletable, match="only major version 3 is supported"):
        delete_latest_upload(graph["uploads"][1])


def test_metadata_disagreeing_with_the_database_is_refused(make_graph):
    """Two independent derivations; if they differ, one assumption is wrong - refuse."""
    graph = make_graph()
    # Declare an observation upload 2 did not create.
    _write_metadata(
        graph["target_dir"] / "upload_2",
        version_number=2,
        observations=[("Peel-x0001", 1), ("Peel-x0001", 2), ("Peel-x0002", 2)],
        canon_sites=[("cs1", 1)],
        conformer_sites=[("csc1", 1)],  # the v2 conf site: created, but not declared
        xtalform_sites=[("xs1", 1)],
    )

    with pytest.raises(NotDeletable, match="disagree about what upload 2 created"):
        delete_latest_upload(graph["uploads"][1])

    # Nothing was touched.
    assert SiteObservation.objects.filter(pk=graph["newcomer_obs"].pk).exists()


def test_unmatched_declaration_warns_rather_than_failing(make_graph):
    """The loader skips entries it cannot process, so the bundle may over-declare."""
    graph = make_graph()
    # Delete the row but leave the declaration: the bundle now claims more than exists.
    graph["newcomer_obs"].delete()
    Experiment.objects.filter(pk=graph["newcomer"].pk).delete()

    plan = delete_latest_upload(graph["uploads"][1])

    assert any("matched no row" in w for w in plan.warnings)
    assert not SiteObservation.objects.filter(version=2).exists()


# --------------------------------------------------------------------------- #
# The deletion itself
# --------------------------------------------------------------------------- #


def test_delete_latest_upload_restores_the_previous_state(make_graph):
    """The whole point: after deleting upload 2, the target looks like it did at 1."""
    graph = make_graph()
    v1, v2 = graph["observations"][1], graph["observations"][2]

    plan = delete_latest_upload(graph["uploads"][1])

    # 1. every version-2 row is gone, every version-1 row survives
    assert not SiteObservation.objects.filter(pk=v2.pk).exists()
    assert not SiteObservation.objects.filter(pk=graph["newcomer_obs"].pk).exists()
    assert SiteObservation.objects.filter(pk=v1.pk).exists()

    # 2. the superseded flag is recomputed from what remains
    v1.refresh_from_db()
    assert v1.superseded is False

    # 3. the pose survives, its main is the version-1 row, its name is untouched
    pose = Pose.objects.get(pk=graph["pose"].pk)
    assert pose.main_site_observation_id == v1.pk
    assert pose.display_name == "Peel-x0001_A_501"

    # 4. the pose whose only member was a version-2 row is gone
    assert not Pose.objects.filter(pk=graph["doomed_pose"].pk).exists()

    # 5. the surviving CanonSiteConf kept its row, with the dangling reference cleared -
    #    the CanonSite -> XtalformSite -> v1 cascade never fired
    conf_site = CanonSiteConf.objects.get(pk=graph["conf_site"].pk)
    assert conf_site.ref_site_observation_id is None
    assert CanonSite.objects.filter(pk=graph["canon_site"].pk).exists()
    assert XtalformSite.objects.filter(pk=graph["xtalform_site"].pk).exists()

    # 6. the experiment created by upload 2 is gone; the carried-over one survives,
    #    re-pointed at upload 1
    assert not Experiment.objects.filter(pk=graph["newcomer"].pk).exists()
    carried = Experiment.objects.get(pk=graph["carried"].pk)
    assert carried.experiment_upload_id == graph["uploads"][0].pk

    # 7. the upload row itself is gone
    assert not ExperimentUpload.objects.filter(pk=graph["uploads"][1].pk).exists()
    assert plan.upload_version == 2


def test_relinked_observation_is_moved_back_to_its_predecessor(make_graph):
    """A survivor re-linked to a site row this upload created must not cascade away.

    The loader recomputes ``canon_site_conf`` on every load, so an observation from an
    earlier upload can end up pointing at a brand-new conformer site. That FK is CASCADE:
    deleting the new row blindly would destroy the survivor.
    """
    graph = make_graph()
    v1 = graph["observations"][1]
    assert v1.canon_site_conf_id == graph["new_conf_site"].pk

    plan = delete_latest_upload(graph["uploads"][1])

    v1.refresh_from_db()
    assert v1.canon_site_conf_id == graph["conf_site"].pk
    assert not CanonSiteConf.objects.filter(pk=graph["new_conf_site"].pk).exists()
    assert (
        "CanonSiteConf",
        v1.pk,
        graph["new_conf_site"].pk,
        graph["conf_site"].pk,
    ) in plan.children_repointed


def test_rows_that_were_already_empty_are_left_alone(make_graph):
    """The orphan sweep removes what the deletion emptied - not what was empty already.

    A target can legitimately carry site rows that no observation points at: A71EV2A has
    10 such CanonSiteConfs from its very first upload. An earlier version of the sweep
    deleted any empty row it found, so peeling an unrelated upload quietly took those with
    it and the target did not come back to its previous state.
    """
    graph = make_graph()
    stranded = CanonSiteConf.objects.create(
        name="never-referenced", canon_site=graph["canon_site"], residues=[], version=1
    )
    lonely_xtalform_site = XtalformSite.objects.create(
        xtalform_site_id="xs-unused",
        xtalform=graph["xtalform_site"].xtalform,
        canon_site=graph["canon_site"],
        lig_chain=CHAIN,
        residues=[],
        version=1,
    )

    plan = delete_latest_upload(graph["uploads"][1])

    assert CanonSiteConf.objects.filter(pk=stranded.pk).exists()
    assert XtalformSite.objects.filter(pk=lonely_xtalform_site.pk).exists()
    assert stranded.pk not in plan.orphaned.get("CanonSiteConf", [])


def test_compounds_the_upload_orphaned_are_swept(make_graph):
    """A compound the deleted upload introduced, now pointed at by nothing, goes.

    Compounds are project-scoped rather than target-scoped, so they are not covered by
    the per-target probes the site models use; without this they survived every deletion
    and accumulated as debris on each load/delete cycle.
    """
    graph = make_graph()
    newcomer_compound = Compound.objects.create(
        smiles="CCN", inchi="InChI=1S/C2H7N", project=graph["project"]
    )
    obs = graph["newcomer_obs"]
    obs.cmpd = newcomer_compound
    obs.save(update_fields=["cmpd"])

    plan = delete_latest_upload(graph["uploads"][1])

    assert not Compound.objects.filter(pk=newcomer_compound.pk).exists()
    assert newcomer_compound.pk in plan.orphaned.get("Compound", [])
    # The compound the surviving observation still uses is untouched.
    assert Compound.objects.filter(pk=graph["compound"].pk).exists()


def test_a_compound_still_referenced_elsewhere_is_kept(make_graph):
    """Anything still pointing at it - here an ExperimentCompound link - saves it."""
    graph = make_graph()
    shared = Compound.objects.create(
        smiles="CCN", inchi="InChI=1S/C2H7N", project=graph["project"]
    )
    obs = graph["newcomer_obs"]
    obs.cmpd = shared
    obs.save(update_fields=["cmpd"])
    # A surviving experiment also claims it.
    ExperimentCompound.objects.create(experiment=graph["carried"], compound=shared)

    plan = delete_latest_upload(graph["uploads"][1])

    assert Compound.objects.filter(pk=shared.pk).exists()
    assert shared.pk not in plan.orphaned.get("Compound", [])


def test_a_compound_already_orphaned_is_left_alone(make_graph):
    """Only what THIS deletion orphaned - the rule the CanonSiteConf bug taught."""
    graph = make_graph()
    stranded = Compound.objects.create(
        smiles="CCCO", inchi="InChI=1S/C3H8O", project=graph["project"]
    )

    delete_latest_upload(graph["uploads"][1])

    assert Compound.objects.filter(pk=stranded.pk).exists()


def test_the_deleted_upload_directory_is_removed_entirely(make_graph):
    """Nothing of the deleted upload may survive on disk - directory included.

    The carried-over experiment's pdb_info was rewritten into upload 2 when that upload
    re-supplied the crystal. Deleting upload 2 restores the path from upload 1's
    metadata, which frees the file, which lets the whole directory go.
    """
    graph = make_graph()
    target_dir = graph["target_dir"]
    upload_1, upload_2 = target_dir / "upload_1", target_dir / "upload_2"
    carried = graph["carried"]

    overwritten = upload_2 / "crystallographic_files" / "Peel-x0001" / "Peel-x0001.pdb"
    debris = upload_2 / "unreferenced.map"
    assert overwritten.is_file() and debris.is_file()
    assert str(carried.pdb_info).endswith(
        "upload_2/crystallographic_files/Peel-x0001/Peel-x0001.pdb"
    )

    plan = delete_latest_upload(graph["uploads"][1])

    # The experiment now points at upload 1's copy, which is still on disk.
    carried.refresh_from_db()
    assert "upload_1/crystallographic_files/Peel-x0001/Peel-x0001.pdb" in str(
        carried.pdb_info
    )
    assert (graph["media_root"] / str(carried.pdb_info)).is_file()
    assert (carried.pk, ["pdb_info", "pdb_info_source_file"]) in (
        plan.experiment_files_restored
    )

    # And with nothing referencing it, the entire directory is gone.
    assert not upload_2.exists()

    # Upload 1 is untouched, and its metadata - the input the next deletion needs - is
    # never swept even though no database row references it.
    assert (upload_1 / "Peel-x0001_v1.pdb").is_file()
    assert (upload_1 / "meta_aligner.yaml").is_file()

    # The bundle archive goes. It sits in the target directory, which is NOT where
    # ExperimentUpload.file resolves to, so the obvious lookup finds nothing - that is
    # how six real bundles totalling 8.7 GB survived a full peel.
    assert not (target_dir / "bundle_upload_2.tgz").is_file()
    # ...while the surviving upload's archive is untouched.
    assert (target_dir / "bundle_upload_1.tgz").is_file()


def test_a_file_this_upload_first_supplied_is_cleared(make_graph):
    """When no surviving upload declares the file, the previous value was *nothing*.

    SoakDB registers a crystal as an Experiment long before XCA aligns it, so an
    experiment can predate the upload that first supplies its crystallographic files -
    real case: A71EV2A-p0411, absent from uploads 1-5 and first supplied by upload 6.
    Clearing the field restores the earlier state exactly, and frees the file so the
    directory can go. Leaving it would strand the row on a directory being removed.
    """
    graph = make_graph()
    target_dir = graph["target_dir"]
    # Take upload 1's declaration away, so upload 2 is the first to supply the file.
    meta = yaml.safe_load((target_dir / "upload_1" / "meta_aligner.yaml").read_text())
    del meta["crystals"]["Peel-x0001"]["crystallographic_files"]
    (target_dir / "upload_1" / "meta_aligner.yaml").write_text(yaml.safe_dump(meta))

    plan = delete_latest_upload(graph["uploads"][1])

    carried = graph["carried"]
    carried.refresh_from_db()
    assert not carried.pdb_info
    assert (carried.pk, ["pdb_info", "pdb_info_source_file"]) in (
        plan.experiment_files_cleared
    )
    # and with nothing referencing it any more, the directory still goes
    assert not (target_dir / "upload_2").exists()


def test_peeling_leaves_no_deleted_upload_directories(make_graph):
    """After peeling back to the first upload, only its directory remains on disk."""
    graph = make_graph(uploads=3)
    target = graph["target"]
    target_dir = graph["target_dir"]

    for _ in range(2):
        upload = (
            ExperimentUpload.objects.filter(target=target)
            .order_by("upload_version")
            .last()
        )
        delete_latest_upload(upload)

    assert (target_dir / "upload_1").is_dir()
    assert not (target_dir / "upload_2").exists()
    assert not (target_dir / "upload_3").exists()
    # and the surviving experiment points into the directory that is left
    graph["carried"].refresh_from_db()
    assert "upload_1/" in str(graph["carried"].pdb_info)
    assert (graph["media_root"] / str(graph["carried"].pdb_info)).is_file()


def test_empty_upload_directory_is_removed(make_graph):
    """With nothing left referencing it, the directory itself goes."""
    graph = make_graph()
    # Point the carried-over experiment back at its own upload's file, so upload 2 holds
    # nothing a survivor needs.
    carried = graph["carried"]
    carried.pdb_info = "target_loader_data/Peel_proposal/upload_1/Peel-x0001_v1.pdb"
    carried.save(update_fields=["pdb_info"])

    delete_latest_upload(graph["uploads"][1])

    assert not (graph["target_dir"] / "upload_2").exists()
    assert (graph["target_dir"] / "upload_1").is_dir()


# --------------------------------------------------------------------------- #
# Eligibility
# --------------------------------------------------------------------------- #


def test_only_the_latest_upload_can_be_deleted(make_graph):
    graph = make_graph(uploads=3)

    with pytest.raises(NotDeletable, match="Only the latest upload"):
        delete_latest_upload(graph["uploads"][0])


def test_the_only_upload_cannot_be_deleted(make_graph):
    """There is no previous state to restore to - delete_target is the right tool."""
    graph = make_graph(uploads=1)

    with pytest.raises(NotDeletable, match="deleting the target itself"):
        delete_latest_upload(graph["uploads"][0])


# --------------------------------------------------------------------------- #
# Iteration
# --------------------------------------------------------------------------- #


def test_deleting_repeatedly_peels_uploads_off_one_at_a_time(make_graph):
    """Deleting upload N makes N-1 the latest, and the invariants hold at every step."""
    graph = make_graph(uploads=3)
    target = graph["target"]

    for expected_version in (3, 2):
        upload = (
            ExperimentUpload.objects.filter(target=target)
            .order_by("upload_version")
            .last()
        )
        assert upload.upload_version == expected_version
        delete_latest_upload(upload)

        remaining = ExperimentUpload.objects.filter(target=target)
        assert remaining.count() == expected_version - 1

        # No experiment points at an upload that is gone.
        assert (
            not Experiment.objects.filter(experiment_upload__target=target)
            .exclude(experiment_upload__in=remaining)
            .exists()
        )

        # Every surviving pose has a main that is one of its own members and is current.
        for pose in Pose.objects.filter(canon_site=graph["canon_site"]):
            main = pose.main_site_observation
            assert main is not None
            assert main.pose_id == pose.pk
            assert main.superseded is False

        # The metadata of the uploads that remain is still on disk for the next pass.
        for remaining_upload in remaining:
            metadata = (
                graph["target_dir"]
                / remaining_upload.upload_data_dir
                / "meta_aligner.yaml"
            )
            assert metadata.is_file()

    # Down to the first upload, which is not deletable.
    last = ExperimentUpload.objects.get(target=target)
    assert last.upload_version == 1
    assert (
        SiteObservation.objects.filter(experiment__experiment_upload=last).count() == 1
    )
    with pytest.raises(NotDeletable, match="deleting the target itself"):
        delete_latest_upload(last)


def test_experiments_are_repointed_to_their_earliest_surviving_observation(make_graph):
    """``min(surviving version)``, not "the previous upload" - the rule that iterates.

    A pre-fix target can have an experiment stamped with the latest upload while holding
    observations from much earlier ones. Re-pointing it one step back would merely delete
    it on the next pass; re-pointing it at the upload matching its earliest surviving
    observation is a fixed point.
    """
    graph = make_graph(uploads=3)
    carried = graph["carried"]
    # Mis-stamp the experiment the way the loader used to, wholesale, before 6314c86b.
    carried.experiment_upload = graph["uploads"][2]
    carried.save(update_fields=["experiment_upload"])

    delete_latest_upload(graph["uploads"][2])

    carried.refresh_from_db()
    # Its earliest surviving observation is version 1, so it belongs to upload 1 - not to
    # upload 2, which is where a naive "previous upload" rule would have put it.
    assert carried.experiment_upload_id == graph["uploads"][0].pk


# --------------------------------------------------------------------------- #
# Superseding
# --------------------------------------------------------------------------- #


def test_restore_superseded_is_a_pure_function_of_surviving_rows(make_graph):
    """No undo journal: the flag is derivable from whichever rows remain."""
    graph = make_graph(uploads=3)
    target = graph["target"]

    # Scramble the flags; they must come back from the data alone.
    SiteObservation.objects.filter(experiment__experiment_upload__target=target).update(
        superseded=True
    )
    changed = restore_superseded(target)

    assert changed["SiteObservation"] > 0
    flags = {
        obs.version: obs.superseded
        for obs in SiteObservation.objects.filter(experiment=graph["carried"])
    }
    assert flags == {1: True, 2: True, 3: False}


# --------------------------------------------------------------------------- #
# Dry run
# --------------------------------------------------------------------------- #


def test_plan_upload_deletion_changes_nothing(make_graph):
    graph = make_graph()
    before = {
        "observations": SiteObservation.objects.count(),
        "poses": Pose.objects.count(),
        "uploads": ExperimentUpload.objects.count(),
        "main": Pose.objects.get(pk=graph["pose"].pk).main_site_observation_id,
    }

    plan = plan_upload_deletion(graph["uploads"][1])

    assert sorted(plan.doomed["SiteObservation"]) == sorted(
        [graph["observations"][2].pk, graph["newcomer_obs"].pk]
    )
    assert plan.poses_repointed == [
        (graph["pose"].pk, graph["observations"][2].pk, graph["observations"][1].pk)
    ]
    assert SiteObservation.objects.count() == before["observations"]
    assert Pose.objects.count() == before["poses"]
    assert ExperimentUpload.objects.count() == before["uploads"]
    assert (
        Pose.objects.get(pk=graph["pose"].pk).main_site_observation_id == before["main"]
    )
    assert (graph["target_dir"] / "upload_2" / "unreferenced.map").is_file()


# --------------------------------------------------------------------------- #
# Management command
# --------------------------------------------------------------------------- #


def test_command_deletes_the_latest_upload_when_asked(make_graph, capsys):
    graph = make_graph()

    call_command("delete_target", "--title", "Peel", "--latest-upload")

    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 1
    assert "Deleted upload 2 of Peel" in capsys.readouterr().out


def test_bare_invocation_prints_help_and_deletes_nothing(make_graph, capsys):
    """A destructive command must do nothing, and explain itself, on a bare run."""
    graph = make_graph()

    call_command("delete_target")

    out = capsys.readouterr().out
    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 2
    assert "--latest-upload" in out
    assert "--entire-target" in out
    # the help, not a listing of the database
    assert "Uploads:" not in out


def test_naming_a_target_without_an_action_lists_its_uploads(
    make_graph, make_project, make_target, capsys
):
    graph = make_graph()
    other = make_target(make_project("other-proposal"), title="Elsewhere")

    call_command("delete_target", "--pk", str(graph["target"].pk))

    out = capsys.readouterr().out
    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 2
    assert "Uploads:" in out
    assert "upload 1" in out and "upload 2" in out
    assert f"--pk {graph['target'].pk} --latest-upload" in out
    # scoped to the target that was named
    assert other.title not in out


def test_title_alone_lists_every_target_with_that_title(
    make_graph, make_project, make_target, capsys
):
    """Listing is not destructive, so an ambiguous title is answered, not refused.

    Showing each match with its pk and proposal is how a caller finds the one they
    want; the delete paths still refuse rather than guess.
    """
    graph = make_graph()
    twin = make_target(make_project("other-proposal"), title="Peel")

    call_command("delete_target", "--title", "Peel")

    out = capsys.readouterr().out
    assert "2 targets match" in out
    assert f"--pk {graph['target'].pk}" in out
    assert f"--pk {twin.pk}" in out
    assert "proposal other-proposal" in out
    assert "--pk to pick one" in out
    # still nothing deleted
    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 2


def test_title_alone_with_one_match_lists_just_it(make_graph, capsys):
    graph = make_graph()

    call_command("delete_target", "--title", "Peel")

    out = capsys.readouterr().out
    assert "targets match" not in out
    assert f"--pk {graph['target'].pk} --latest-upload" in out


def test_listing_an_unknown_title_says_so(db):
    with pytest.raises(CommandError, match="No target titled 'Nope'"):
        call_command("delete_target", "--title", "Nope")


def test_an_action_without_a_target_is_refused(make_graph):
    make_graph()

    with pytest.raises(CommandError, match="Name the target to act on"):
        call_command("delete_target", "--latest-upload")


def test_command_refuses_two_actions_at_once(make_graph):
    make_graph()

    with pytest.raises(CommandError, match="different operations"):
        call_command(
            "delete_target", "--title", "Peel", "--latest-upload", "--entire-target"
        )


def test_command_refuses_pk_and_title_together(make_graph):
    graph = make_graph()

    with pytest.raises(CommandError, match="either --pk or --title"):
        call_command(
            "delete_target",
            "--pk",
            str(graph["target"].pk),
            "--title",
            "Peel",
            "--latest-upload",
        )


def test_command_dry_run_reports_without_changing_anything(make_graph, capsys):
    graph = make_graph()

    call_command("delete_target", "--title", "Peel", "--latest-upload", "--dry-run")

    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 2
    out = capsys.readouterr().out
    assert "Dry run - nothing was changed." in out
    assert "delete" in out


def test_command_entire_target_removes_everything(make_graph):
    graph = make_graph()

    call_command("delete_target", "--title", "Peel", "--entire-target")

    assert not ExperimentUpload.objects.filter(target=graph["target"]).exists()
    assert not SiteObservation.objects.filter(experiment__isnull=False).exists()


def test_command_refuses_in_production(make_graph, settings):
    make_graph()
    settings.DEPLOYMENT_MODE = "PRODUCTION"

    with pytest.raises(CommandError, match="disabled in production"):
        call_command("delete_target", "--title", "Peel", "--latest-upload")


def test_command_reports_an_unknown_target(db):
    with pytest.raises(CommandError, match="No target titled 'Nope'"):
        call_command("delete_target", "--title", "Nope", "--latest-upload")


def test_command_refuses_an_ambiguous_title(make_graph, make_project, make_target):
    """A title is unique only within a proposal, so a bare title can name several.

    The command is destructive, so it must not pick one: it lists the candidates and
    makes the caller narrow it down.
    """
    make_graph()
    other = make_project("other-proposal")
    twin = make_target(other, title="Peel")

    with pytest.raises(CommandError, match="names 2 targets") as excinfo:
        call_command("delete_target", "--title", "Peel", "--latest-upload")

    message = str(excinfo.value)
    assert "--proposal" in message
    assert f"--pk {twin.pk}" in message
    assert "proposal other-proposal" in message


def test_command_disambiguates_by_proposal(make_graph, make_project, make_target):
    """--proposal picks the right one out of identically titled targets."""
    graph = make_graph()
    make_target(make_project("other-proposal"), title="Peel")

    call_command(
        "delete_target", "--title", "Peel", "--proposal", "proposal", "--latest-upload"
    )

    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 1


def test_command_accepts_a_primary_key(make_graph, make_project, make_target):
    """--pk addresses the target directly - the form that never needs --proposal."""
    graph = make_graph()
    make_target(make_project("other-proposal"), title="Peel")

    call_command("delete_target", "--pk", str(graph["target"].pk), "--latest-upload")

    assert ExperimentUpload.objects.filter(target=graph["target"]).count() == 1


def test_command_reports_an_unknown_primary_key(db):
    with pytest.raises(CommandError, match="No target with --pk 999"):
        call_command("delete_target", "--pk", "999", "--latest-upload")


def test_command_refuses_a_single_upload_target(make_graph):
    make_graph(uploads=1)

    with pytest.raises(CommandError, match="deleting the target itself") as excinfo:
        call_command("delete_target", "--title", "Peel", "--latest-upload")

    # This is the one refusal --entire-target actually answers, so it is offered.
    assert "--entire-target" in str(excinfo.value)


def test_other_refusals_do_not_suggest_entire_target(make_graph):
    """Offering --entire-target for, say, unreadable metadata would be wrong advice."""
    graph = make_graph()
    (graph["target_dir"] / "upload_2" / "meta_aligner.yaml").unlink()

    with pytest.raises(CommandError, match="cannot establish") as excinfo:
        call_command("delete_target", "--title", "Peel", "--latest-upload")

    assert "--entire-target" not in str(excinfo.value)
