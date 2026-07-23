"""Tests for the small lookup helpers in ``hypothesis.definitions`` (#991).

``IntTypes`` and ``VectTypes`` are pure in-memory tables mapping external names
to two-letter codes. The tests live under ``viewer/tests/`` because that is the
suite's single ``testpaths`` root (as ``test_api_utils.py`` already does).
"""
import pytest

from hypothesis.definitions import IntTypes, VectTypes


def test_get_int_conv_known_version_and_name():
    int_types = IntTypes()
    assert int_types.get_int_conv("DE", "hbond") == "HB"
    assert int_types.get_int_conv("PR", "vdW") == "VD"


def test_get_int_conv_unknown_version_returns_none():
    int_types = IntTypes()
    assert int_types.get_int_conv("ZZ", "hbond") is None


def test_get_int_conv_unknown_name_raises():
    int_types = IntTypes()
    with pytest.raises(KeyError):
        int_types.get_int_conv("DE", "nonsense")


def test_define_int_types_returns_defaults():
    int_types = IntTypes()
    ver_choices, default_ver, type_choices, default_type = int_types.define_int_types()
    assert default_ver == "DE"
    assert default_type == "UK"
    assert ("DE", "Default") in ver_choices
    assert ("HB", "H-bond") in type_choices


@pytest.mark.parametrize(
    "name,code",
    [
        ("additions", "AD"),
        ("deletions", "DE"),
        ("linkers", "LI"),
        ("ring", "RI"),
    ],
)
def test_translate_vect_types(name, code):
    assert VectTypes().translate_vect_types(name) == code


def test_translate_vect_types_unknown_raises():
    with pytest.raises(KeyError):
        VectTypes().translate_vect_types("unknown")


def test_get_vect_types_returns_table():
    assert ("AD", "Addition") in VectTypes().get_vect_types()
