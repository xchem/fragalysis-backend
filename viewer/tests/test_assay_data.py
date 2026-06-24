"""Tests for the pure cell parsers in ``viewer.assay_data`` (issue #991).

``process_float_value`` / ``process_int_value`` parse a single spreadsheet cell
into a ``(raw, modifier, value, parsing_error, error_msg)`` 5-tuple, recognising
the ``<``, ``>``, ``<=``, ``>=`` modifier prefixes and flagging unparseable
cells. ``process_text_value`` and ``get_unit`` are likewise pure. The DataFrame
/ model-backed wrappers (``process_float`` etc.) need ResultValue* rows and are
deferred per the issue.
"""
import math

import pytest

from viewer.assay_data import (
    get_unit,
    process_float_value,
    process_int_value,
    process_text_value,
)


@pytest.mark.parametrize(
    "value,modifier,parsed",
    [
        ("5.0", None, "5.0"),
        ("<=5.0", "<=", "5.0"),
        (">=5.0", ">=", "5.0"),
        ("<5.0", "<", "5.0"),
        (">5.0", ">", "5.0"),
        (".5", None, ".5"),
        ("-3.2", None, "-3.2"),
    ],
)
def test_process_float_value_parses(value, modifier, parsed):
    raw, mod, val, parsing_error, error_msg = process_float_value(value)
    assert raw == value
    assert mod == modifier
    assert val == parsed
    assert parsing_error is False
    assert math.isnan(error_msg)


def test_process_float_value_unparsable():
    raw, mod, val, parsing_error, error_msg = process_float_value("abc")
    assert raw == "abc"
    assert mod is None
    assert val is None
    assert parsing_error is True
    assert error_msg == "Unable to parse abc to float"


@pytest.mark.parametrize(
    "value,modifier,parsed",
    [
        ("5", None, "5"),
        ("<=5", "<=", "5"),
        (">=5", ">=", "5"),
        ("<5", "<", "5"),
        (">5", ">", "5"),
        ("-3", None, "-3"),
    ],
)
def test_process_int_value_parses(value, modifier, parsed):
    raw, mod, val, parsing_error, error_msg = process_int_value(value)
    assert raw == value
    assert mod == modifier
    assert val == parsed
    assert parsing_error is False
    assert math.isnan(error_msg)


def test_process_int_value_unparsable():
    raw, mod, val, parsing_error, error_msg = process_int_value("abc")
    assert raw == "abc"
    assert mod is None
    assert val is None
    assert parsing_error is True
    assert error_msg == "Unable to parse abc to int"


def test_process_text_value():
    assert process_text_value("hello") == ("hello", "hello", None)


@pytest.mark.parametrize(
    "title,expected",
    [
        ("IC50 (nM)", "nM"),
        ("pIC50", ""),
        ("Activity (uM) extra", "uM"),
    ],
)
def test_get_unit(title, expected):
    assert get_unit(title) == expected
