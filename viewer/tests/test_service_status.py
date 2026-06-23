"""Regression tests for service_status probe error handling (issue #980).

When neo4j is unreachable the ``fragmentation_graph`` probe used to let the
neo4j ``ServiceUnavailable`` exception propagate out of the probe. APScheduler
then logged a full, multi-frame traceback every 30s (see the stack-0 logs after
the h/w outage of 18 June 2026). The probe must instead catch that exception,
log a single concise line, and report ``DEGRADED``.
"""
import logging
from contextlib import contextmanager

from neo4j.exceptions import ServiceUnavailable

from service_status import services
from service_status.utils import State


class _FakeSession:
    """A neo4j session whose ``run`` always raises the given exception."""

    def __init__(self, exc):
        self._exc = exc

    def run(self, *_args, **_kwargs):
        raise self._exc


class _FakeDriver:
    """A neo4j driver handing out sessions that fail on ``run``."""

    def __init__(self, exc):
        self._exc = exc

    @contextmanager
    def session(self):
        yield _FakeSession(self._exc)


# The probe is wrapped by @service_query (which touches the DB); exercise the
# underlying function directly via the reference functools.wraps leaves behind.
fragmentation_graph = services.fragmentation_graph.__wrapped__

_UNAVAILABLE = ServiceUnavailable("Couldn't connect to graph.graph-x.svc:7687")


def test_fragmentation_graph_degraded_when_graph_unavailable(monkeypatch):
    """An unreachable graph is reported as DEGRADED, not raised."""
    monkeypatch.setattr(
        services, "get_driver", lambda **_kwargs: _FakeDriver(_UNAVAILABLE)
    )

    # Must not raise, and must report DEGRADED.
    assert fragmentation_graph() == State.DEGRADED


def test_fragmentation_graph_logs_concise_message(monkeypatch, caplog):
    """The unreachable graph is logged as one concise line, no traceback."""
    monkeypatch.setattr(
        services, "get_driver", lambda **_kwargs: _FakeDriver(_UNAVAILABLE)
    )

    # The 'service_status' logger is configured with propagate=False, so its
    # records never reach caplog's root handler. Enable propagation for the
    # duration of the test (monkeypatch restores it) so caplog can see them.
    service_logger = logging.getLogger("service_status")
    monkeypatch.setattr(service_logger, "propagate", True)

    with caplog.at_level("WARNING", logger="service_status"):
        fragmentation_graph()

    # A single record naming the failure (matching the issue's desired text),
    # and crucially no exception traceback attached to it.
    messages = [record.getMessage() for record in caplog.records]
    assert any("ServiceUnavailable" in message for message in messages)
    assert all(record.exc_info is None for record in caplog.records)
