# infections.py
#
# A utility that provides a list of infections (forced internal errors).
# Infections are injected into the application via the environment variable
# 'INFECTIONS', a comma-separated list of infection names.

from typing import Dict, Set

from django.conf import settings

from api.utils import deployment_mode_is_production

# The built-in set of infections.
# Must be lowercase, but user can use any case in the environment variable.
# Define every name as a constant, and add it and a description to the _CATALOGUE.
INFECTION_STRUCTURE_DOWNLOAD: str = 'structure-download'

# The index is the short-form name of the infection, and the value is the
# description of the infection.
_CATALOGUE: Dict[str, str] = {
    INFECTION_STRUCTURE_DOWNLOAD: 'An error in the DownloadStructures view'
}


def have_infection(name: str) -> bool:
    """Returns True if we've been given the named infection.
    Infections are never present in production mode."""
    return False if deployment_mode_is_production() else name in _get_infections()


def _get_infections() -> Set[str]:
    if settings.INFECTIONS == '':
        return set()
    infections: set[str] = {
        infection
        for infection in settings.INFECTIONS.split(',')
        if infection in _CATALOGUE
    }
    return infections
