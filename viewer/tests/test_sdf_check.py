"""Tests for ``viewer.sdf_check.add_warning`` (#991).

It builds a ``validate_dict`` of parallel ``molecule_name`` / ``field`` /
``warning_string`` lists, which the computed-set upload uses to collect
warnings.
"""
from viewer import sdf_check


def empty_validate_dict():
    return {"molecule_name": [], "field": [], "warning_string": []}


def test_add_warning_appends_parallel_lists():
    vd = empty_validate_dict()
    result = sdf_check.add_warning("mol1", "field1", "boom", vd)
    assert result["molecule_name"] == ["mol1"]
    assert result["field"] == ["field1"]
    assert result["warning_string"] == ["boom"]
