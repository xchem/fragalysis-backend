"""Tests for the computed-set download endpoint after the ComputedMolecule ->
SiteObservation unification.

``GET /api/compound-sets/{pk}/download/`` (``ComputedSetView.download``) zips up
each computed observation's user-uploaded pdb. Post-unification those pdbs live on
``SiteObservation.virtual_pdb_info`` and only some observations have one - the
rest reference an existing experimental observation instead. This exercises:

- an observation *with* an uploaded pdb is included, named via
  ``SiteObservation.get_filename()`` (the auto-assigned ``#<hash>`` suffix stripped);
- an observation *without* one is skipped entirely - no file and, crucially, no
  spurious ``<name>_MISSING`` entry (the in-loop ``virtual_pdb_info`` guard).

The test settings leave ``TA_AUTH_SERVICE`` unset, so proposal membership is
resolved from ``Project.user_id`` and ``DEPLOYMENT_MODE`` defaults to DEVELOPMENT.
"""

import zipfile
from io import BytesIO

from viewer.models import ComputedSet, ComputedSetSiteObservation, SiteObservation

# The `db` fixture is requested for its side effect (DB access) only, and
# fixtures legitimately reuse their names as arguments - both standard pytest
# patterns that pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument


def test_download_zips_only_observations_with_uploaded_pdb(
    settings, tmp_path, user, make_project, make_target, authenticated_client
):
    """The download contains exactly the uploaded pdb (correctly named), and the
    referenced observation produces neither a file nor a `_MISSING` entry."""
    media_root = tmp_path / "media"
    (media_root / settings.COMPUTED_SET_MEDIA_DIRECTORY).mkdir(parents=True)
    settings.MEDIA_ROOT = str(media_root)

    project = make_project("proposal", members=[user])
    target = make_target(project, title="DownloadMe")

    # written_sdf_filename left None so the (unscoped) history-SDF branch of the
    # view yields nothing and the zip holds only the pdb we care about here.
    cset = ComputedSet.objects.create(
        name="cset", target=target, owner_user=user, written_sdf_filename=None
    )

    # Observation WITH an uploaded pdb. The stored name carries the auto-assigned
    # '#<hash>' suffix that get_filename() must strip down to "A0486a.pdb".
    pdb_rel = f"{settings.COMPUTED_SET_MEDIA_DIRECTORY}/A0486a#abc.pdb_def"
    (media_root / pdb_rel).write_text("PDBDATA")
    uploaded = SiteObservation.objects.create(smiles="CCO", virtual_pdb_info=pdb_rel)

    # Observation WITHOUT one (it referenced an existing observation on upload).
    referenced = SiteObservation.objects.create(smiles="CCC")

    ComputedSetSiteObservation.objects.create(
        computed_set=cset, site_observation=uploaded
    )
    ComputedSetSiteObservation.objects.create(
        computed_set=cset, site_observation=referenced
    )

    response = authenticated_client.get(f"/api/compound-sets/{cset.pk}/download/")

    assert response.status_code == 200

    with zipfile.ZipFile(BytesIO(response.content)) as archive:
        names = archive.namelist()
        # Exactly the uploaded pdb, under its stripped filename.
        assert names == ["A0486a.pdb"]
        assert archive.read("A0486a.pdb") == b"PDBDATA"
        # The referenced observation must not leak a `_MISSING` placeholder.
        assert not any(name.endswith("_MISSING") for name in names)
