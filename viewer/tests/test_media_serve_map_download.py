"""Tests for map_download's extension->field resolution (issue #984, phase 3).

``media_serve.views.map_download`` picks which SiteObservation file field to
serve by substring-parsing the requested file name (``2fofc``/``fofc``/``event``
=> ``sigmaa_info``/``diff_info``/``event_info``). The parsing is fragile, so the
field selection is pinned here. We patch ``get_response`` to capture the chosen
field rather than hit the DB - the download authorisation itself is covered by
test_security_static_files.py.
"""
# pylint: disable=redefined-outer-name
import pytest
from django.test import RequestFactory

from api.security import ISPyBSafeStaticFiles
from media_serve import views


@pytest.fixture
def captured_field(monkeypatch):
    """Patch get_response to return the field_name map_download selected."""
    monkeypatch.setattr(
        ISPyBSafeStaticFiles, "get_response", lambda self: self.field_name
    )

    def _call(file_path):
        request = RequestFactory().get("/maps/" + file_path)
        return views.map_download(request, file_path)

    return _call


@pytest.mark.parametrize(
    "file_path,expected_field",
    [
        ("Mpro-x0001_2fofc.map", "sigmaa_info"),
        ("Mpro-x0001_fofc.map", "diff_info"),
        ("Mpro-x0001_event.map", "event_info"),
    ],
)
def test_map_download_selects_field_by_extension(
    captured_field, file_path, expected_field
):
    assert captured_field(file_path) == expected_field


def test_map_download_unknown_extension_selects_no_field(captured_field):
    """A name with no recognised map extension resolves to no field (None)."""
    assert captured_field("Mpro-x0001_apo.map") is None
