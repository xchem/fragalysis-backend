"""Tests for the compound-set SDF validators in ``viewer.sdf_check`` (#991).

These validators build a ``validate_dict`` of parallel ``molecule_name`` /
``field`` / ``warning_string`` lists, appending a record for every problem
found. The tests construct RDKit mols from SMILES (with properties set via
``SetProp``) - no database or filesystem is involved.
"""
import pytest
from rdkit import Chem

from viewer import sdf_check


def empty_validate_dict():
    return {"molecule_name": [], "field": [], "warning_string": []}


def make_mol(name="mol1", props=None):
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", name)
    for key, value in (props or {}).items():
        mol.SetProp(key, value)
    return mol


def warnings_count(validate_dict):
    return len(validate_dict["warning_string"])


def test_add_warning_appends_parallel_lists():
    vd = empty_validate_dict()
    result = sdf_check.add_warning("mol1", "field1", "boom", vd)
    assert result["molecule_name"] == ["mol1"]
    assert result["field"] == ["field1"]
    assert result["warning_string"] == ["boom"]


def test_check_sdf_illegal_filename_warns():
    vd = sdf_check.check_sdf("compound-set_thing.txt", empty_validate_dict())
    assert warnings_count(vd) == 1
    assert vd["field"] == ["_File_name"]


def test_check_sdf_valid_filename_no_warning():
    vd = sdf_check.check_sdf("compound-set_thing.sdf", empty_validate_dict())
    assert warnings_count(vd) == 0


@pytest.mark.parametrize("name", ["good-name_1.0", "abc.def"])
def test_check_name_characters_legal(name):
    vd = sdf_check.check_name_characters(name, empty_validate_dict())
    assert warnings_count(vd) == 0


def test_check_name_characters_illegal():
    vd = sdf_check.check_name_characters("bad name!", empty_validate_dict())
    # space and '!' are both illegal -> two warnings
    assert warnings_count(vd) == 2


def test_check_url_valid_returns_none():
    assert sdf_check.check_url("https://example.com") is None


def test_check_url_invalid_returns_false():
    assert sdf_check.check_url("not a url") is False


def test_missing_field_check_present():
    mol = make_mol(props={"ref_mols": "x0107"})
    vd = sdf_check.missing_field_check(mol, "ref_mols", empty_validate_dict())
    assert warnings_count(vd) == 0


def test_missing_field_check_absent():
    mol = make_mol()
    vd = sdf_check.missing_field_check(mol, "ref_mols", empty_validate_dict())
    assert warnings_count(vd) == 1
    assert vd["field"] == ["ref_mols"]


def test_check_blank_prop_warns_on_empty_and_bad_url():
    mol = make_mol(props={"submitter_name": "", "ref_url": "not a url"})
    vd = sdf_check.check_blank_prop(mol, empty_validate_dict())
    # empty submitter_name -> 1 warning; bad ref_url -> 1 warning
    assert warnings_count(vd) == 2


def test_check_blank_prop_ignores_listed_props():
    mol = make_mol(props={"ref_pdb": ""})
    vd = sdf_check.check_blank_prop(mol, empty_validate_dict())
    assert warnings_count(vd) == 0


def test_check_mol_props_missing_ref_pdb_and_lhs_pdb():
    mol = make_mol(props={"ref_mols": "x0107"})
    vd = sdf_check.check_mol_props(mol, empty_validate_dict())
    # ref_mols present, but neither ref_pdb nor lhs_pdb -> 1 warning
    assert warnings_count(vd) == 1
    assert vd["field"] == ["ref_pdb/lhs_pdb"]


def test_check_mol_props_ref_pdb_satisfies():
    mol = make_mol(props={"ref_mols": "x0107", "ref_pdb": "x0107.pdb"})
    vd = sdf_check.check_mol_props(mol, empty_validate_dict())
    assert warnings_count(vd) == 0


def test_check_mol_props_lhs_pdb_satisfies():
    mol = make_mol(props={"ref_mols": "x0107", "lhs_pdb": "x0107.pdb"})
    vd = sdf_check.check_mol_props(mol, empty_validate_dict())
    assert warnings_count(vd) == 0


def test_check_ver_name_matches():
    mol = make_mol(name="ver_1.2")
    vd = sdf_check.check_ver_name(mol, "ver_1.2", empty_validate_dict())
    assert warnings_count(vd) == 0


def test_check_ver_name_mismatch():
    mol = make_mol(name="ver_1.0")
    vd = sdf_check.check_ver_name(mol, "ver_1.2", empty_validate_dict())
    assert warnings_count(vd) == 1


def test_check_compound_set_no_generation_date():
    mol = make_mol()
    vd = sdf_check.check_compound_set(mol, empty_validate_dict())
    assert warnings_count(vd) == 1
    assert "no generation_date" in vd["warning_string"][0]


def test_check_compound_set_malformed_generation_date():
    mol = make_mol(props={"generation_date": "2020/01/01"})
    vd = sdf_check.check_compound_set(mol, empty_validate_dict())
    assert warnings_count(vd) == 1


def test_check_compound_set_valid_generation_date():
    mol = make_mol(props={"generation_date": "2020-01-01"})
    vd = sdf_check.check_compound_set(mol, empty_validate_dict())
    assert warnings_count(vd) == 0
