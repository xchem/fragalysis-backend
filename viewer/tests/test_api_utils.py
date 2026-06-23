"""Tests for the small cross-cutting helpers in ``api.utils`` (issue #984, phase 3).

``validate_tas`` is the Target-Access-String gate used wherever a TAS is
accepted, and ``deployment_mode_is_production`` switches strict behaviour on -
both are high-traffic and previously untested.

Note: ``api.utils.TAS_REGEX_RE`` is compiled from ``settings.TAS_REGEX`` at
import time, so these tests use the default pattern (``^(lb|sw)\\d{5}-\\d+$``)
rather than overriding the setting (which would not recompile the regex).
"""
import pytest

from api.utils import deployment_mode_is_production, validate_tas


@pytest.mark.parametrize("tas", ["lb12345-1", "sw00000-99", "lb54321-123"])
def test_validate_tas_accepts_well_formed(tas):
    valid, error = validate_tas(tas)
    assert valid is True
    assert error is None


@pytest.mark.parametrize(
    "tas",
    [
        "",  # empty
        "xx12345-1",  # wrong prefix
        "lb1234-1",  # too few digits
        "lb12345",  # missing visit
        "lb12345-",  # missing visit number
        " lb12345-1",  # leading space
    ],
)
def test_validate_tas_rejects_malformed(tas, settings):
    valid, error = validate_tas(tas)
    assert valid is False
    assert error == settings.TAS_REGEX_ERROR_MSG


def test_deployment_mode_is_production_true(settings):
    settings.DEPLOYMENT_MODE = "PRODUCTION"
    assert deployment_mode_is_production() is True


def test_deployment_mode_is_production_false(settings):
    settings.DEPLOYMENT_MODE = "DEVELOPMENT"
    assert deployment_mode_is_production() is False
