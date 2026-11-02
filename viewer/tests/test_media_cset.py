# pylint: disable=redefined-outer-name,unused-argument
"""Tests for computed-set (RHS) uploads through the media storage seam (issue 1038).

A computed-set upload is a chain of Celery tasks that can run on different
workers, so the files pass between them as storage names, never local paths:

    UploadComputedSetView -> save_tmp_file()   tmp/cset/<uuid>/<file>
    validate_compound_set                     reads the SDF from storage
    process_compound_set -> PdbOps            stores each PDB in the zip
                         -> MolOps            moves the PDBs, stores the
                                              virtual .mol files, moves the SDF

Each test runs with both media storages (the ``media_backend`` fixture).
"""
import importlib
import io
import re
import zipfile
from datetime import datetime
from datetime import timezone as dt_timezone
from pathlib import Path

import pytest
from django.apps import apps
from django.core.files.uploadedfile import SimpleUploadedFile
from rdkit import Chem

from viewer import media_storage
from viewer.cset_upload import MolOps, PdbOps
from viewer.models import ComputedSet, ExperimentUpload
from viewer.tasks import validate_compound_set
from viewer.utils import save_tmp_file

_TEST_SDF = Path(__file__).resolve().parents[2] / "tests/test_data/compund-set_test.sdf"
_MOL_BLOCK = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))


def _zip_bytes(files: dict[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    return buffer.getvalue()


@pytest.fixture
def loaded_target(user, make_project, make_target):
    """A target with one experiment upload (MolOps needs its directory)."""
    target = make_target(make_project("proposal-cset", members=[user]), title="Mpro")
    target.zip_archive = "Mpro_lb1"
    target.save()
    ExperimentUpload.objects.create(
        project=target.project,
        target=target,
        file="Mpro_lb1.tgz",
        upload_data_dir="upload_1",
        commit_datetime=datetime(2026, 1, 1, tzinfo=dt_timezone.utc),
        committer=user,
    )
    return target


def _mol_ops(target, sdf_name="tmp/cset/u1/my set.sdf", zfile=None):
    return MolOps(
        user_id=None,
        sdf_filename=sdf_name,
        submitter_name="",
        submitter_method="",
        target=target.pk,
        version="ver_1.2",
        zfile=zfile or {},
        zfile_hashvals={},
        computed_set_id=None,
    )


# The upload ----------------------------------------------------------------


def test_save_tmp_file_stores_the_upload_and_returns_its_name(media_backend):
    name = save_tmp_file(SimpleUploadedFile("my set.sdf", b"sdf"))

    assert re.match(r"^tmp/cset/[0-9a-f]{32}/my_set.sdf$", name)
    assert media_storage.read_bytes(name) == b"sdf"


# validate_compound_set ------------------------------------------------------


def test_validate_reports_a_missing_sdf(media_backend):
    _, validate_dict, validated, params = validate_compound_set(
        {"user_id": None, "sdf_file": "tmp/cset/u1/missing.sdf", "target": 1}
    )

    assert not validated
    assert validate_dict["molecule_name"] == []
    assert params["sdf"] == "tmp/cset/u1/missing.sdf"


def test_validate_reads_the_stored_sdf(media_backend, mocker):
    name = media_storage.store_bytes("tmp/cset/u1/set.sdf", _TEST_SDF.read_bytes())
    check = mocker.patch(
        "viewer.tasks.check_compound_set", side_effect=lambda mol, d, **_: d
    )

    _, _, _, params = validate_compound_set(
        {"user_id": None, "sdf_file": name, "target": 1}
    )

    blank_mol = check.call_args.args[0]
    assert blank_mol.GetProp("submitter_name") == "Alessandro Contini"
    assert params["sdf"] == name


# PdbOps ----------------------------------------------------------------------


@pytest.mark.django_db
def test_pdb_zip_files_are_stored_and_named(media_backend):
    zip_name = media_storage.store_bytes(
        "tmp/cset/u1/refs.zip",
        _zip_bytes({"refs/Mpro-x0001.pdb": b"ATOM 1\n", "refs/readme.txt": b"x"}),
    )

    zfile, hashvals = PdbOps().save_pdb_zip(zip_name)

    assert zfile is not None
    assert hashvals == {}
    assert list(zfile) == ["Mpro-x0001"]
    assert re.match(r"^tmp/cset/[0-9a-f]{32}/Mpro-x0001.pdb$", zfile["Mpro-x0001"])
    assert media_storage.read_bytes(zfile["Mpro-x0001"]) == b"ATOM 1\n"


def test_no_pdb_zip_gives_no_files(media_backend):
    assert PdbOps().save_pdb_zip(None) == (None, None)


# MolOps: PDBs ---------------------------------------------------------------


def test_process_pdb_moves_the_pdb_to_the_computed_set_directory(
    media_backend, loaded_target
):
    tmp_pdb = media_storage.store_bytes("tmp/cset/u2/Mpro-x0001.pdb", b"ATOM\n")
    zfile = {"Mpro-x0001": tmp_pdb}
    ops = _mol_ops(loaded_target, zfile=zfile)

    name = ops.process_pdb("Mpro-x0001", zfile, {})

    assert re.match(r"^computed_set_data/Mpro-x0001.pdb_[0-9a-f]{32}$", name)
    assert media_storage.read_bytes(name) == b"ATOM\n"
    assert not media_storage.exists(tmp_pdb)


def test_process_pdb_gives_one_name_for_a_pdb_used_twice(media_backend, loaded_target):
    """Two molecules can reference the same PDB. Both get the moved file."""
    zfile = {
        "Mpro-x0001": media_storage.store_bytes("tmp/cset/u2/Mpro-x0001.pdb", b"A")
    }
    ops = _mol_ops(loaded_target, zfile=zfile)

    first = ops.process_pdb("Mpro-x0001", zfile, {})
    second = ops.process_pdb("Mpro-x0001", zfile, {})

    assert first == second
    assert media_storage.read_bytes(second) == b"A"


def test_process_pdb_of_an_unknown_code_is_none(media_backend, loaded_target):
    assert _mol_ops(loaded_target).process_pdb("Mpro-x9999", {}, {}) is None


def test_delete_uploaded_pdbs_removes_stored_and_moved_pdbs(
    media_backend, loaded_target
):
    zfile = {
        "a": media_storage.store_bytes("tmp/cset/u2/a.pdb", b"a"),
        "b": media_storage.store_bytes("tmp/cset/u2/b.pdb", b"b"),
    }
    ops = _mol_ops(loaded_target, zfile=zfile)
    moved = ops.process_pdb("a", zfile, {})

    ops.delete_uploaded_pdbs()

    assert not media_storage.exists(moved)
    assert not media_storage.exists(zfile["b"])


# MolOps: virtual .mol files ---------------------------------------------------


def test_virtual_mol_files_are_stored_together(media_backend, loaded_target):
    ops = _mol_ops(loaded_target)
    name = str(ops.virtual_root / "cset_upload_1_v1a_A0486a_1.mol")

    ops.write_virtual_mol(name, _MOL_BLOCK)
    assert not media_storage.exists(name)  # Not stored until the end

    ops.store_virtual_files()

    assert name.startswith("target_loader_data/Mpro_lb1/upload_1/virtual_files/")
    assert media_storage.read_text(name) == _MOL_BLOCK


def test_read_virtual_mol_finds_new_and_stored_files(media_backend, loaded_target):
    ops = _mol_ops(loaded_target)
    stored = media_storage.store_bytes("computed_set_data/old.mol", _MOL_BLOCK)
    new = str(ops.virtual_root / "new.mol")
    ops.write_virtual_mol(new, _MOL_BLOCK)

    for name in (stored, new):
        mol = ops.read_virtual_mol(name)
        assert Chem.MolToSmiles(mol) == "CCO"


# MolOps: the SDF -------------------------------------------------------------


def test_store_sdf_moves_the_sdf_and_records_storage_names(
    media_backend, loaded_target, user
):
    sdf = media_storage.store_bytes("tmp/cset/u1/my_set.sdf", b"sdf")
    computed_set = ComputedSet.objects.create(
        name="cset", md_ordinal=2, target=loaded_target, owner_user=user
    )

    _mol_ops(loaded_target, sdf_name=sdf).store_sdf(computed_set)

    computed_set.refresh_from_db()
    assert computed_set.submitted_sdf == "my_set.sdf"
    assert computed_set.written_sdf_filename == (
        "computed_set_data/cset_upload_2_my_set.sdf"
    )
    assert media_storage.read_bytes(computed_set.written_sdf_filename) == b"sdf"
    assert not media_storage.exists(sdf)


# The download -----------------------------------------------------------------


@pytest.mark.parametrize("legacy_absolute", [False, True])
def test_download_includes_the_written_sdf(
    media_backend, settings, user, loaded_target, authenticated_client, legacy_absolute
):
    name = media_storage.store_bytes("computed_set_data/cset_upload_1_a.sdf", b"SDF")
    stored_value = str(Path(settings.MEDIA_ROOT) / name) if legacy_absolute else name
    cset = ComputedSet.objects.create(
        name="cset",
        target=loaded_target,
        owner_user=user,
        submitted_sdf="a.sdf",
        written_sdf_filename=stored_value,
    )

    response = authenticated_client.get(f"/api/compound-sets/{cset.pk}/download/")

    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        assert archive.read("a.sdf") == b"SDF"


# The migration -----------------------------------------------------------------


def test_migration_strips_the_media_root_from_written_sdf_filenames(
    settings, user, make_project, make_target
):
    migration = importlib.import_module(
        "viewer.migrations.0174_written_sdf_filename_storage_relative"
    )
    settings.MEDIA_ROOT = "/code/media/"
    target = make_target(make_project("proposal-mig", members=[user]))
    absolute = ComputedSet.objects.create(
        name="a",
        target=target,
        owner_user=user,
        written_sdf_filename="/code/media/computed_set_data/a.sdf",
    )
    relative = ComputedSet.objects.create(
        name="b",
        target=target,
        owner_user=user,
        written_sdf_filename="computed_set_data/b.sdf",
    )
    unset = ComputedSet.objects.create(name="c", target=target, owner_user=user)

    migration.to_storage_names(apps, None)

    for cset in (absolute, relative, unset):
        cset.refresh_from_db()
    assert absolute.written_sdf_filename == "computed_set_data/a.sdf"
    assert relative.written_sdf_filename == "computed_set_data/b.sdf"
    assert unset.written_sdf_filename is None
    history = ComputedSet.history.filter(name="a").values_list(
        "written_sdf_filename", flat=True
    )
    assert set(history) == {"computed_set_data/a.sdf"}

    migration.to_absolute_paths(apps, None)

    absolute.refresh_from_db()
    assert absolute.written_sdf_filename == "/code/media/computed_set_data/a.sdf"
    assert set(
        ComputedSet.history.filter(name="b").values_list(
            "written_sdf_filename", flat=True
        )
    ) == {"/code/media/computed_set_data/b.sdf"}
