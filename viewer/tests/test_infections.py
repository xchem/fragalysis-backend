"""Tests for the infection (forced-error) catalogue (issue #984, phase 3).

``api.infections`` parses the ``INFECTIONS`` env-driven setting into a set of
recognised infection names and gates them behind the deployment mode -
infections must never fire in production. These are pure, DB-free helpers.
"""
from api import infections
from api.infections import INFECTION_STRUCTURE_DOWNLOAD, have_infection


def test_empty_setting_yields_no_infections(settings):
    settings.INFECTIONS = ""
    assert infections._get_infections() == set()  # pylint: disable=protected-access


def test_known_infections_are_parsed(settings):
    settings.INFECTIONS = INFECTION_STRUCTURE_DOWNLOAD
    assert infections._get_infections() == {  # pylint: disable=protected-access
        INFECTION_STRUCTURE_DOWNLOAD
    }


def test_unknown_infection_is_filtered_out(settings):
    settings.INFECTIONS = f"{INFECTION_STRUCTURE_DOWNLOAD},not-a-real-infection"
    assert infections._get_infections() == {  # pylint: disable=protected-access
        INFECTION_STRUCTURE_DOWNLOAD
    }


def test_have_infection_true_when_present_and_not_production(settings):
    settings.DEPLOYMENT_MODE = "DEVELOPMENT"
    settings.INFECTIONS = INFECTION_STRUCTURE_DOWNLOAD
    assert have_infection(INFECTION_STRUCTURE_DOWNLOAD) is True


def test_have_infection_always_false_in_production(settings):
    """Even if requested, an infection never fires in production mode."""
    settings.DEPLOYMENT_MODE = "PRODUCTION"
    settings.INFECTIONS = INFECTION_STRUCTURE_DOWNLOAD
    assert have_infection(INFECTION_STRUCTURE_DOWNLOAD) is False
