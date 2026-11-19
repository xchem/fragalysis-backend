"""Regression tests for connecting to a TLS-only graph (issue #1039).

The graph now requires an encrypted connection, reached by hostname using a
``bolt+s://`` URI (e.g. ``bolt+s://graph-x.xchem-dev.diamond.ac.uk:7687``).
``NEO4J_QUERY`` must therefore accept a complete URI, which has to reach the
neo4j driver unchanged, while a bare hostname keeps the original plain
``bolt://<host>:7687`` connection.

These tests exercise the real ``fragutils`` driver construction, replacing only
the neo4j ``GraphDatabase.driver`` call so no graph is contacted.
"""
from service_status import services
from service_status.utils import State

# The probe is wrapped by @service_query (which touches the DB); exercise the
# underlying function directly via the reference functools.wraps leaves behind.
fragmentation_graph = services.fragmentation_graph.__wrapped__


def test_secure_graph_uri_reaches_driver_unchanged(settings, driver_calls):
    settings.NEO4J_QUERY = "bolt+s://graph-x.xchem-dev.diamond.ac.uk:7687"
    settings.NEO4J_AUTH = "neo4j/secret"

    assert fragmentation_graph() == State.OK

    assert driver_calls == [
        (
            ("bolt+s://graph-x.xchem-dev.diamond.ac.uk:7687",),
            {"auth": ("neo4j", "secret")},
        )
    ]


def test_bare_graph_hostname_uses_plain_bolt(settings, driver_calls):
    settings.NEO4J_QUERY = "graph.graph-a.svc"
    settings.NEO4J_AUTH = "neo4j/secret"

    assert fragmentation_graph() == State.OK

    assert driver_calls == [
        (("bolt://graph.graph-a.svc:7687",), {"auth": ("neo4j", "secret")})
    ]
