# pylint: disable=redefined-outer-name,unused-argument
"""Tests for the models and serializers that read media files (issue 1038).

They read through ``viewer.media_storage``, so each test runs with both media
storages (the ``media_backend`` fixture). The models are not saved: these
methods use only the instances' fields.
"""
from pathlib import Path

import pytest

from viewer import media_storage
from viewer.models import Experiment, ExperimentUpload, SiteObservation, Target
from viewer.serializers import (
    ProtMapInfoSerializer,
    ProtPDBBoundInfoSerializer,
    ProtPDBInfoSerializer,
    TargetSerializer,
)
from viewer.target_loader import XTALFORMS_FILE

_PREFIX = "target_loader_data/Mpro_lb1/upload_1"
_ASSEMBLIES = """\
assemblies:
  dimer:
    reference: Mpro-x0001
    biomol: A,B
"""
_REF_PDB = Path(
    "target_loader_data/Mpro_lb1/upload_1/crystallographic_files/Mpro-x0001/Mpro-x0001.pdb"
)


@pytest.fixture
def upload():
    return ExperimentUpload(
        target=Target(zip_archive="Mpro_lb1"), upload_data_dir="upload_1"
    )


# ExperimentUpload ------------------------------------------------------------


def test_upload_prefix_is_the_storage_name_of_the_upload_directory(upload):
    assert upload.upload_prefix() == _PREFIX


# SiteObservation -------------------------------------------------------------


def test_get_ligand_mol_file_reads_the_stored_file(media_backend):
    media_storage.store_bytes("target_loader_data/a.mol", "mol block\n")

    assert SiteObservation(
        ligand_mol="target_loader_data/a.mol"
    ).get_ligand_mol_file() == ("mol block\n")


@pytest.mark.parametrize("ligand_mol", [None, "", "target_loader_data/missing.mol"])
def test_get_ligand_mol_file_is_empty_without_a_file(media_backend, ligand_mol):
    assert SiteObservation(ligand_mol=ligand_mol).get_ligand_mol_file() == ""


def test_virtual_ligand_mol_name_is_the_storage_name():
    obs = SiteObservation(virtual_ligand_mol="computed_set_data/A0486a.mol")

    assert obs.virtual_ligand_mol_name == "computed_set_data/A0486a.mol"


@pytest.mark.parametrize("value", [None, ""])
def test_virtual_ligand_mol_name_is_none_when_unset(value):
    assert SiteObservation(virtual_ligand_mol=value).virtual_ligand_mol_name is None


# Template protein (TargetSerializer) ----------------------------------------


def test_template_protein_path_is_the_reference_pdb(media_backend, upload):
    media_storage.store_bytes(f"{_PREFIX}/{XTALFORMS_FILE}", _ASSEMBLIES)
    media_storage.store_bytes(str(_REF_PDB), "ATOM\n")

    assert TargetSerializer().get_template_protein_path(upload) == _REF_PDB


def test_template_protein_path_is_none_without_the_reference_pdb(media_backend, upload):
    media_storage.store_bytes(f"{_PREFIX}/{XTALFORMS_FILE}", _ASSEMBLIES)

    assert TargetSerializer().get_template_protein_path(upload) is None


def test_template_protein_path_is_none_without_the_assemblies_file(
    media_backend, upload
):
    assert TargetSerializer().get_template_protein_path(upload) is None


@pytest.mark.parametrize(
    "assemblies", ["other: {}\n", "assemblies: {}\n", "assemblies:\n  dimer: {}\n"]
)
def test_template_protein_path_is_none_for_a_bad_assemblies_file(
    media_backend, upload, assemblies
):
    media_storage.store_bytes(f"{_PREFIX}/{XTALFORMS_FILE}", assemblies)

    assert TargetSerializer().get_template_protein_path(upload) is None


# File content serializers ---------------------------------------------------


def test_pdb_data_is_the_stored_pdb(media_backend):
    media_storage.store_bytes("target_loader_data/a.pdb", "ATOM\n")
    obs = SiteObservation(experiment=Experiment(pdb_info="target_loader_data/a.pdb"))

    assert ProtPDBInfoSerializer().get_pdb_data(obs) == "ATOM\n"


def test_bound_pdb_data_is_the_stored_file(media_backend):
    media_storage.store_bytes("target_loader_data/bound.pdb", "HETATM\n")

    obs = SiteObservation(bound_file="target_loader_data/bound.pdb")

    assert ProtPDBBoundInfoSerializer().get_bound_pdb_data(obs) == "HETATM\n"
    assert ProtPDBBoundInfoSerializer().get_bound_pdb_data(SiteObservation()) is None


def test_map_data_is_the_stored_file(media_backend):
    media_storage.store_bytes("target_loader_data/event.ccp4", "map\n")

    obs = SiteObservation(event_file="target_loader_data/event.ccp4")

    assert ProtMapInfoSerializer().get_map_data(obs) == "map\n"
    assert ProtMapInfoSerializer().get_map_data(SiteObservation()) is None
