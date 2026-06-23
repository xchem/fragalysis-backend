"""Tests for service_status.utils.services enable/disable/toggle (issue #984, phase 3).

``services()`` is the manual override that flips a Service row between
NOT_CONFIGURED (disabled) and DEGRADED (enabled, pending first probe). The
enable/disable/toggle branches (including a name in *both* lists) had no tests.
ServiceState rows are seeded by the service_status migrations.
"""
# pylint: disable=redefined-outer-name
import pytest

from service_status.models import Service, ServiceState
from service_status.utils import State, services


def _make_service(name: str, state: State) -> Service:
    return Service.objects.create(
        service=name,
        display_name=name,
        last_state=ServiceState.objects.get(state=state),
    )


@pytest.mark.django_db
def test_enable_promotes_not_configured_to_degraded():
    _make_service("ispyb", State.NOT_CONFIGURED)

    services(enable=["ispyb"])

    assert Service.objects.get(service="ispyb").last_state_id == State.DEGRADED


@pytest.mark.django_db
def test_enable_leaves_already_active_service_untouched():
    """Enabling an OK service must not clobber it back to DEGRADED."""
    _make_service("ispyb", State.OK)

    services(enable=["ispyb"])

    assert Service.objects.get(service="ispyb").last_state_id == State.OK


@pytest.mark.django_db
def test_disable_sets_not_configured():
    _make_service("ispyb", State.OK)

    services(disable=["ispyb"])

    assert Service.objects.get(service="ispyb").last_state_id == State.NOT_CONFIGURED


@pytest.mark.django_db
def test_toggle_when_in_both_lists_enables_a_disabled_service():
    """A name in both enable and disable toggles from its current state."""
    _make_service("ispyb", State.NOT_CONFIGURED)

    services(enable=["ispyb"], disable=["ispyb"])

    assert Service.objects.get(service="ispyb").last_state_id == State.DEGRADED


@pytest.mark.django_db
def test_toggle_when_in_both_lists_disables_an_enabled_service():
    _make_service("ispyb", State.OK)

    services(enable=["ispyb"], disable=["ispyb"])

    assert Service.objects.get(service="ispyb").last_state_id == State.NOT_CONFIGURED


@pytest.mark.django_db
def test_unknown_service_is_skipped_without_error():
    """An unknown service name is logged and skipped, not raised."""
    # No Service rows exist; this must not raise.
    services(enable=["does-not-exist"], disable=["also-missing"])

    assert Service.objects.count() == 0


@pytest.mark.django_db
def test_none_arguments_are_treated_as_empty():
    _make_service("ispyb", State.OK)

    services(enable=None, disable=None)

    # Nothing requested, so state is unchanged.
    assert Service.objects.get(service="ispyb").last_state_id == State.OK
