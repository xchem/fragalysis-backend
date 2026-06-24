"""Tests for service_status.init_services activation logic (issue #982).

`init_services()` reconciles the Service table against the
`ENABLE_SERVICE_STATUS` setting. It used to only *demote* services missing from
the list to NOT_CONFIGURED; it never *promoted* a requested service out of
NOT_CONFIGURED. Because a Service row defaults to NOT_CONFIGURED and the probe
self-skips in that state, a requested service could stay NOT_CONFIGURED forever.
These tests pin down that a requested service is activated, while an already
active one is left untouched, and an unrequested one is disabled.
"""
# Fixtures legitimately reuse their own name as a test argument, are requested
# for their side effect only, and reach for a module "constant" deliberately -
# all standard pytest patterns that pylint misreads (cf. conftest.py).
# pylint: disable=redefined-outer-name,unused-argument,protected-access
import pytest

from service_status import utils
from service_status.models import Service, ServiceState
from service_status.utils import State, init_services


def _set_state(service_name: str, state: State, display_name: str = "x") -> None:
    """Create (or reset) a Service row in a known state."""
    Service.objects.update_or_create(
        service=service_name,
        defaults={
            "display_name": display_name,
            "last_state": ServiceState.objects.get(state=state),
        },
    )


@pytest.fixture
def as_service_check_host(monkeypatch, settings):
    """Run init_services as the check host, without a real scheduler."""
    monkeypatch.setattr(utils, "_HOSTNAME", utils._SERVICE_CHECK_HOSTNAME)
    settings.SERVICE_STATUS_SCHEDULER_ENABLED = False
    settings.ENABLE_SERVICE_STATUS = "fragmentation_graph:ispyb:keycloak:squonk"


@pytest.mark.django_db
def test_requested_not_configured_service_is_activated(as_service_check_host):
    """A requested service stuck at NOT_CONFIGURED is promoted to DEGRADED."""
    _set_state("fragmentation_graph", State.NOT_CONFIGURED, "Fragmentation graph")

    init_services()

    service = Service.objects.get(service="fragmentation_graph")
    assert service.last_state_id == State.DEGRADED


@pytest.mark.django_db
def test_already_active_requested_service_is_preserved(as_service_check_host):
    """An active requested service is not reset (e.g. OK must stay OK)."""
    _set_state("ispyb", State.OK, "Access control (ISPyB)")

    init_services()

    service = Service.objects.get(service="ispyb")
    assert service.last_state_id == State.OK


@pytest.mark.django_db
def test_unrequested_service_is_disabled(as_service_check_host):
    """A service absent from ENABLE_SERVICE_STATUS is set NOT_CONFIGURED."""
    # 'discourse' is defined in services.py but not in the requested list.
    _set_state("discourse", State.OK, "Discourse")

    init_services()

    service = Service.objects.get(service="discourse")
    assert service.last_state_id == State.NOT_CONFIGURED
