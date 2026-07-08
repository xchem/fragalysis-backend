"""Tests for the target-deletion feature.

Covers the reusable service function (``viewer.target_delete.delete_target``) and
the guarded ``DELETE /api/targets/{pk}/`` endpoint on ``TargetView``:

- the service removes a target's DB graph and media files, for that target alone;
- the endpoint is rejected in production, for anonymous users and for non-members,
  and succeeds (204) for a member on a DEVELOPMENT instance.

The test settings leave ``TA_AUTH_SERVICE`` unset, so proposal membership is
resolved from ``Project.user_id`` (no external service involved), and
``DEPLOYMENT_MODE`` defaults to ``DEVELOPMENT``.
"""

from pathlib import Path

from viewer.models import (
    ComputedSet,
    ComputedSetSiteObservation,
    SiteObservation,
    Target,
)
from viewer.target_delete import delete_target

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


def _add_computed_set(target, user, media_root, settings, name="cset") -> ComputedSet:
    """Create a ComputedSet for ``target`` with on-disk files, returning it.

    Lays down the three kinds of computed-set file the deleter must remove
    individually: the written SDF (absolute path in a TextField), the submitted
    SDF (FileField under computed_set_data/) and a computed observation's uploaded
    pdb file (SiteObservation.virtual_pdb_info).
    """
    cset_dir = media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY

    written_sdf = cset_dir / f"{name}_written.sdf"
    written_sdf.write_text("written")
    submitted_rel = f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/{name}_submitted.sdf"
    (media_root / submitted_rel).write_text("submitted")

    computed_set = ComputedSet.objects.create(
        name=name,
        target=target,
        submitted_sdf=submitted_rel,
        written_sdf_filename=str(written_sdf),
        spec_version=1.0,
        owner_user=user,
    )

    pdb_rel = f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/{name}_mol.pdb"
    (media_root / pdb_rel).write_text("pdb")
    # A computed observation is now a SiteObservation (its FKs are all nullable,
    # so no experiment graph is needed here), linked to the set via the
    # ComputedSetSiteObservation join table.
    obs = SiteObservation.objects.create(virtual_pdb_info=pdb_rel)
    ComputedSetSiteObservation.objects.create(
        computed_set=computed_set, site_observation=obs
    )
    return computed_set


# --------------------------------------------------------------------------- #
# Service function
# --------------------------------------------------------------------------- #


def test_delete_target_removes_db_rows_and_media(
    db, settings, tmp_path, user, make_project, make_target
):
    """delete_target removes the target, its computed sets/molecules and its
    media (loader subdir + computed-set files), leaving other targets alone."""
    media_root = _make_media_root(settings, tmp_path)
    loader_dir = media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY

    project = make_project("proposal", members=[user])

    # Target under test, with a loader subdirectory and a computed set.
    target = make_target(project, title="DeleteMe")
    target.zip_archive = "DeleteMe_proposal"
    target.save()
    target_subdir = loader_dir / "DeleteMe_proposal"
    target_subdir.mkdir()
    (target_subdir / "data.pdb").write_text("x")
    cset = _add_computed_set(target, user, media_root, settings, name="goner")

    # An unrelated target whose data must survive untouched.
    other = make_target(project, title="KeepMe")
    other.zip_archive = "KeepMe_proposal"
    other.save()
    other_subdir = loader_dir / "KeepMe_proposal"
    other_subdir.mkdir()
    (other_subdir / "data.pdb").write_text("y")
    other_cset = _add_computed_set(other, user, media_root, settings, name="keeper")

    delete_target(target)

    # Target and its computed graph are gone.
    assert not Target.objects.filter(pk=target.pk).exists()
    assert not ComputedSet.objects.filter(pk=cset.pk).exists()
    # The computed-set/observation join rows go with the ComputedSet cascade.
    assert not ComputedSetSiteObservation.objects.filter(computed_set=cset).exists()

    # Its media is gone: loader subdir and all three computed-set files.
    assert not target_subdir.exists()
    assert not (
        media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY / "goner_written.sdf"
    ).exists()
    assert not (
        media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY / "goner_submitted.sdf"
    ).exists()
    assert not (
        media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY / "goner_mol.pdb"
    ).exists()

    # The shared computed_set_data directory itself is preserved.
    assert (media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY).is_dir()

    # The unrelated target is completely untouched.
    assert Target.objects.filter(pk=other.pk).exists()
    assert ComputedSet.objects.filter(pk=other_cset.pk).exists()
    assert other_subdir.exists()
    assert (
        media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY / "keeper_written.sdf"
    ).exists()
    assert (
        media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY / "keeper_mol.pdb"
    ).exists()


def test_delete_target_without_media_is_safe(
    db, settings, tmp_path, user, make_project, make_target
):
    """A target with no zip_archive and no computed sets deletes cleanly and
    never touches the parent loader directory."""
    media_root = _make_media_root(settings, tmp_path)
    loader_dir = media_root / settings.TARGET_LOADER_MEDIA_DIRECTORY

    target = make_target(make_project("proposal", members=[user]), title="Bare")

    delete_target(target)

    assert not Target.objects.filter(pk=target.pk).exists()
    # A blank zip_archive must not cause the parent dir to be removed.
    assert loader_dir.is_dir()


# --------------------------------------------------------------------------- #
# Endpoint guards
# --------------------------------------------------------------------------- #


def test_delete_endpoint_member_succeeds(
    authenticated_client, user, make_project, make_target
):
    """A member deleting on a DEVELOPMENT instance gets 204 and the target goes."""
    target = make_target(make_project("proposal", members=[user]), title="DeleteMe")

    response = authenticated_client.delete(f"/api/targets/{target.pk}/")

    assert response.status_code == 204
    assert not Target.objects.filter(pk=target.pk).exists()


def test_delete_endpoint_rejected_in_production(
    settings, authenticated_client, user, make_project, make_target
):
    """Deletion is forbidden when the instance is in production mode."""
    settings.DEPLOYMENT_MODE = "PRODUCTION"
    target = make_target(make_project("proposal", members=[user]), title="DeleteMe")

    response = authenticated_client.delete(f"/api/targets/{target.pk}/")

    assert response.status_code == 403
    assert Target.objects.filter(pk=target.pk).exists()


def test_delete_endpoint_rejects_anonymous(api_client, make_project, make_target):
    """An unauthenticated request cannot delete a target."""
    target = make_target(make_project("proposal"), title="DeleteMe")

    response = api_client.delete(f"/api/targets/{target.pk}/")

    assert response.status_code in (401, 403)
    assert Target.objects.filter(pk=target.pk).exists()


def test_delete_endpoint_rejects_non_member(
    authenticated_client, make_project, make_target
):
    """An authenticated non-member of the target's project cannot delete it."""
    # The authenticated user is not a member of this project.
    target = make_target(make_project("members-only"), title="DeleteMe")

    response = authenticated_client.delete(f"/api/targets/{target.pk}/")

    assert response.status_code == 403
    assert Target.objects.filter(pk=target.pk).exists()


def test_delete_endpoint_skips_membership_when_auth_off(
    settings, api_client, make_project, make_target
):
    """With AUTHENTICATE_UPLOAD off, an anonymous request may delete (dev only)."""
    settings.AUTHENTICATE_UPLOAD = False
    target = make_target(make_project("members-only"), title="DeleteMe")

    response = api_client.delete(f"/api/targets/{target.pk}/")

    assert response.status_code == 204
    assert not Target.objects.filter(pk=target.pk).exists()
