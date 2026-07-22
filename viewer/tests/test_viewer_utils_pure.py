"""Pure-function tests for the small helpers in ``viewer.utils`` (issue #991).

These complement ``test_utils.py`` (which covers the filesystem/SDF helpers)
by exercising the string/identifier helpers, the ``alphanumerator`` generator,
``flatten_dict`` and the filesystem-only ``calculate_sha256`` /
``sanitize_directory_name``. None of these need a database; the few that touch
disk use ``tmp_path``.
"""
import itertools

import pytest

from viewer.utils import (
    alphanumerator,
    calculate_sha256,
    clean_object_id,
    flatten_dict,
    longcode_from_tag,
    sanitize_directory_name,
    strip_exp_code,
    strip_version,
)


def test_alphanumerator_from_empty_start():
    """With no start the generator yields a, b, c ... from the beginning."""
    gen = alphanumerator()
    assert list(itertools.islice(gen, 3)) == ["a", "b", "c"]


def test_alphanumerator_drop_first_default():
    """Default drop_first=True returns the value *after* start_from."""
    gen = alphanumerator(start_from="a")
    assert list(itertools.islice(gen, 3)) == ["b", "c", "d"]


def test_alphanumerator_keep_first():
    """drop_first=False makes start_from itself the first value yielded."""
    gen = alphanumerator(start_from="a", drop_first=False)
    assert list(itertools.islice(gen, 3)) == ["a", "b", "c"]


def test_alphanumerator_is_case_insensitive():
    """start_from is lower-cased before matching."""
    gen = alphanumerator(start_from="B", drop_first=False)
    assert next(gen) == "b"


def test_alphanumerator_rolls_over_to_two_letters():
    """After 'z' the generator continues with two-letter words."""
    gen = alphanumerator(start_from="z")
    assert next(gen) == "aa"


@pytest.mark.parametrize(
    "name,expected",
    [
        # The '-x' anchored form: only the tail after '-x' is rewritten.
        ("Mpro-x0107+A+1+1", "Mpro-x0107/A/1/1"),
        ("Mpro-x0107_A_1_1", "Mpro-x0107/A/1/1"),
        # No '-x' anchor: '+' and '_' are rewritten across the whole string.
        ("foo+bar_baz", "foo/bar/baz"),
        # Nothing to replace.
        ("plain", "plain"),
    ],
)
def test_clean_object_id(name, expected):
    assert clean_object_id(name) == expected


@pytest.mark.parametrize(
    "value,expected",
    [
        ("XX01ZVNS2B-x0673/B/501/1", ("XX01ZVNS2B-x0673/B/501", 1)),
        ("a/b/12", ("a/b", 12)),
    ],
)
def test_strip_version(value, expected):
    assert strip_version(value) == expected


def test_strip_version_custom_separator():
    assert strip_version("a-b-7", separator="-") == ("a-b", 7)


def test_longcode_from_tag():
    assert longcode_from_tag("XX-x0673/B/501/1") == "XX-x0673_B_501_v1"


def test_longcode_from_tag_custom_separator():
    assert longcode_from_tag("a-b-2", separator="-") == "a_b_v2"


def test_strip_exp_code():
    """The chunk after the first '-<letter>' marker is returned."""
    assert strip_exp_code("Mpro-x0107") == "0107"


def test_strip_exp_code_non_standard_raises():
    """A code without a '-<letter>' marker raises a ValueError."""
    with pytest.raises(ValueError, match="Non-standard experiment code"):
        strip_exp_code("plaincode")


def test_flatten_dict_shallow():
    """A flat dict (values are not dicts-of-dicts) is yielded unchanged."""
    result = dict(flatten_dict({"a": {"x": 1}, "b": {"y": 2}}))
    assert result == {"a": {"x": 1}, "b": {"y": 2}}


def test_flatten_dict_nested_to_depth():
    """Nested dicts are flattened: a key path of length n yields an
    (n+1)-tuple ending in the leaf value."""
    nested = {"a": {"b": {"x": 1}, "c": {"y": 2}}}
    result = list(flatten_dict(nested, depth=2))
    assert result == [("a", "b", {"x": 1}), ("a", "c", {"y": 2})]


def test_flatten_dict_skips_non_dict_values():
    """Values that are not dicts (no .values()) are skipped, not raised."""
    assert dict(flatten_dict({"a": 1, "b": "text"})) == {}


def test_calculate_sha256(tmp_path):
    """The digest matches the known SHA-256 of the file's bytes."""
    import hashlib

    target = tmp_path / "blob.bin"
    payload = b"fragalysis" * 1000
    target.write_bytes(payload)
    assert calculate_sha256(target) == hashlib.sha256(payload).hexdigest()


def test_sanitize_directory_name_replaces_illegal_chars():
    assert sanitize_directory_name("My Target!/v2") == "My_Target_v2"


def test_sanitize_directory_name_collapses_underscores():
    assert sanitize_directory_name("a   b") == "a_b"


def test_sanitize_directory_name_unique_within_path(tmp_path):
    """A name colliding with an existing dir gets a numeric suffix."""
    (tmp_path / "target").mkdir()
    assert sanitize_directory_name("target", path=tmp_path) == "target_2"
