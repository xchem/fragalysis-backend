"""Unit tests for the pure-Python helpers in ``viewer.utils``.

Ported from the original Django-runner test case to pytest. These need no
database, so they run quickly and in isolation.
"""
import filecmp
import os
from unittest.mock import Mock

import pytest
from django.conf import settings

from viewer.utils import (
    add_prop_to_sdf,
    clean_filename,
    create_media_sub_directory,
    create_squonk_job_request_url,
    delete_media_sub_directory,
    get_https_host,
    is_url,
    word_count,
)


def test_create_squonk_job_request_url():
    assert (
        create_squonk_job_request_url("instance-0000")
        == "data-manager-ui/results/instance/instance-0000"
    )


def test_create_media_sub_directory():
    result = create_media_sub_directory("utils-test-blob")
    try:
        assert result == os.path.join(settings.MEDIA_ROOT, "utils-test-blob")
        assert os.path.isdir(result)
    finally:
        delete_media_sub_directory("utils-test-blob")


def test_delete_media_sub_directory():
    result = create_media_sub_directory("utils-test-blobby")
    assert os.path.isdir(result)
    delete_media_sub_directory("utils-test-blobby")
    assert not os.path.exists(result)


def test_add_prop_to_sdf(tmp_path):
    out_file = str(tmp_path / "out.sdf")
    add_prop_to_sdf(
        "tests/test_data/viewer-utils-test-in.sdf",
        out_file,
        {"TransFSScore": "0.115601"},
    )
    assert os.path.isfile(out_file)
    assert filecmp.cmp("tests/test_data/viewer-utils-test-out.sdf", out_file)


@pytest.mark.parametrize(
    "filepath,expected",
    [
        ("./media/sdfs/Mpro-x3351_0A_rtEVbqf.sdf", "Mpro-x3351_0A.sdf"),
        ("./media/sdfs/Mpro-x3351_0A.sdf", "Mpro-x3351_0A.sdf"),
        ("Mpro-x3351_0A.sdf", "Mpro-x3351_0A.sdf"),
    ],
)
def test_clean_filename(filepath, expected):
    assert clean_filename(filepath) == expected


def test_get_https_host():
    request = Mock()
    request.get_host.return_value = "example.com"
    assert get_https_host(request) == "https://example.com"


@pytest.mark.parametrize(
    "value",
    [
        "https://example.com",
        "http://example.com",
        "ftp://example.com",
        "sftp://example.com",
        "ssh://example.com",
        "file://example.com",
        "ldap://example.com",
    ],
)
def test_is_url_true(value):
    assert is_url(value) is True


@pytest.mark.parametrize("value", ["/data/blob.html", 532, None])
def test_is_url_false(value):
    assert is_url(value) is False


@pytest.mark.parametrize(
    "text,expected",
    [("Hello world", 2), ("Hello", 1), ("", 0), (None, 0)],
)
def test_word_count(text, expected):
    assert word_count(text) == expected
