"""Tests for the tag-file ingestion helpers (issue #984, phase 2).

``viewer.tags.load_tags_from_file`` rewrites tags, poses and quality statuses
from an uploaded metadata file and previously had no tests. The full happy path
needs a completely loaded target (experiments, canon sites, poses) and is best
covered by the loader integration tests; here we pin the pure parsing helpers,
the boolean-column sanitiser (which must flag bad values rather than pass them
silently) and the "not a valid file" guard - the cheap, high-value parts that
guard against silent data corruption.
"""
import pandas as pd
import pytest

from viewer.tags import (
    get_ann_tag,
    get_tag_cols,
    load_tags_from_file,
    sanitize_boolean_column,
    strip_alias,
    strip_catname,
    validate_aliases,
)


@pytest.mark.parametrize(
    "value", ["true", "1", "t", "y", "yes", "True", "TRUE", True, 1]
)
def test_sanitize_boolean_true_representations(value):
    column = pd.Series([value], name="flag")
    result, errors = sanitize_boolean_column(column, [])
    assert result.tolist() == [True]
    assert errors == []


@pytest.mark.parametrize("value", ["false", "0", "n", "no", "False", False, 0, ""])
def test_sanitize_boolean_false_representations(value):
    column = pd.Series([value], name="flag")
    result, errors = sanitize_boolean_column(column, [])
    assert result.tolist() == [False]
    assert errors == []


def test_sanitize_boolean_invalid_value_is_flagged_not_swallowed():
    """An unrecognised value records an error and is passed through unchanged."""
    column = pd.Series(["maybe"], name="flag")
    result, errors = sanitize_boolean_column(column, [])
    assert result.tolist() == ["maybe"]
    assert errors == ["Invalid boolean value 'maybe' in column 'flag'"]


def test_strip_catname_splits_category_and_name():
    assert strip_catname("[Other] P2 Wave 1") == ("Other", "P2 Wave 1")


def test_strip_catname_handles_brackets_in_name():
    assert strip_catname("[Series] foo [bar]") == ("Series", "foo [bar]")


def test_strip_alias_without_prefix_returns_value():
    assert strip_alias("plain") == "plain"


def test_strip_alias_with_prefix_returns_alias_part():
    assert strip_alias("x - alias") == "alias"


def test_get_tag_cols_selects_only_curated_category_columns():
    columns = [
        "Long code",
        "[Series] one",
        "[Forum] two",
        "[Other] three",
        "[CanonSites] auto",  # not a curated category - excluded
    ]
    assert get_tag_cols(columns) == ["[Series] one", "[Forum] two", "[Other] three"]


def test_get_ann_tag_is_stable_md5():
    assert get_ann_tag("Main status") == get_ann_tag("Main status")
    assert get_ann_tag("a") != get_ann_tag("b")


def _alias_frame(rows):
    """Build a DataFrame with the upload-name/alias columns validate_aliases needs."""
    from viewer.tags import TAG_CATEGORIES

    columns: dict[str, list] = {}
    for cat in TAG_CATEGORIES:
        columns[f"{cat} upload name"] = []
        columns[f"{cat} alias"] = []
    df = pd.DataFrame(columns=list(columns.keys()))
    for row in rows:
        df.loc[len(df)] = row
    return df


def test_validate_aliases_accepts_consistent_pairs():
    df = _alias_frame([])
    # One upload name with a single, consistent alias in the first category.
    df.loc[0] = [""] * len(df.columns)
    df.loc[0, "ConformerSites upload name"] = "up1"
    df.loc[0, "ConformerSites alias"] = "alias1"

    result, errors = validate_aliases(df)

    assert errors == []
    assert ("up1", "alias1") in result


def test_validate_aliases_flags_inconsistent_alias():
    df = _alias_frame([])
    for idx in (0, 1):
        df.loc[idx] = [""] * len(df.columns)
        df.loc[idx, "ConformerSites upload name"] = "up1"
    # Same upload name, two different aliases -> inconsistent.
    df.loc[0, "ConformerSites alias"] = "alias1"
    df.loc[1, "ConformerSites alias"] = "alias2"

    _, errors = validate_aliases(df)

    assert any("Inconsistent aliases for tag 'up1'" in e for e in errors)


@pytest.mark.django_db
def test_load_tags_from_file_rejects_non_csv_xlsx(make_project, make_target, tmp_path):
    """A file that is neither CSV nor XLSX yields an error, not a silent pass."""
    target = make_target(make_project("proposal-tags"))
    bad_file = tmp_path / "garbage.bin"
    # Bytes that are invalid UTF-8 (so read_csv raises UnicodeDecodeError) and
    # not a valid workbook (so read_excel raises ValueError).
    bad_file.write_bytes(b"\xff\xfe\x00\x01not a spreadsheet")

    errors = load_tags_from_file(str(bad_file), target)

    assert errors == [f"{bad_file} is not a valid CSV or XLSX file"]
