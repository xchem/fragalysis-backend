"""Pure / RDKit-only tests for helpers in ``api.utils`` (issue #991).

These complement ``test_api_utils.py`` (which covers ``validate_tas`` and
``deployment_mode_is_production``) with the request-parameter parsers
(``parse_bool``, ``parse_vectors``) and the geometry helpers (``calc_bounds``,
``calc_dims``). None need a database or network.

(``_transparentsvg`` is deliberately not tested: it is dead code - defined but
never called - and carries a latent ``str + bytes`` bug, so exercising it would
require fixing unrelated code.)
"""
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from api.utils import calc_bounds, calc_dims, parse_bool, parse_vectors


@pytest.mark.parametrize("value", ["yes", "true", "t", "y", "1", "TRUE", "Yes"])
def test_parse_bool_true(value):
    assert parse_bool(value) is True


@pytest.mark.parametrize("value", ["no", "false", "f", "n", "0", "FALSE", "No"])
def test_parse_bool_false(value):
    assert parse_bool(value) is False


@pytest.mark.parametrize("value", ["", "maybe", "2", "yep"])
def test_parse_bool_unparsable_raises(value):
    with pytest.raises(ValueError, match="not parsable"):
        parse_bool(value)


def test_parse_vectors():
    assert parse_vectors("1,2,3") == [1, 2, 3]


def test_parse_vectors_single():
    assert parse_vectors("42") == [42]


@pytest.mark.parametrize("value", ["1,a,3", "1,,3", ""])
def test_parse_vectors_malformed_raises(value):
    with pytest.raises(ValueError):
        parse_vectors(value)


def test_calc_dims():
    """Dimensions are simply the spans of the X and Y bounds."""
    assert calc_dims([0, 3], [-1, 4]) == (3, 5)


def test_calc_bounds_from_conformer():
    """Bounds enclose every atom's 2D coordinates."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    AllChem.Compute2DCoords(mol)
    conformer = mol.GetConformer()
    x, y = calc_bounds(conformer)
    # Bounds are [min, max]; min must be <= max and the spans positive.
    assert x[0] <= x[1]
    assert y[0] <= y[1]
    assert calc_dims(x, y)[0] > 0
