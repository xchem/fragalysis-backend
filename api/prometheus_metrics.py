"""Prometheus metrics used by the fragalysis API module.
"""
from prometheus_client import Counter


class PrometheusMetrics:
    """A static class to hold the Prometheus metrics for the fragalysis API module.
    Each metric has its own static method to adjust it.
    """

    # Create, and initialise the metrics for this module
    ssh_tunnels = Counter(
        'fragalysis_ssh_tunnels',
        'Number of SSH tunnels successfully created',
    )
    ssh_tunnels.reset()
    ssh_tunnel_failures = Counter(
        'fragalysis_ssh_tunnel_failures',
        'Number of SSH tunnel failures',
    )
    ssh_tunnel_failures.reset()
    ispyb_connections = Counter(
        'fragalysis_ispyb_connections',
        'Number of ISpyB successful connections (excluding retries)',
    )
    ispyb_connections.reset()
    ispyb_connection_attempts = Counter(
        'fragalysis_ispyb_connection_attempts',
        'Number of ISpyB connection retries (after initial failure)',
    )
    ispyb_connection_attempts.reset()
    ispyb_connection_failures = Counter(
        'fragalysis_ispyb_connection_failures',
        'Number of ISpyB connection failures',
    )
    ispyb_connection_failures.reset()
    proposal_cache_hit = Counter(
        'fragalysis_proposal_cache_hit',
        'Number of proposal cache hits',
    )
    proposal_cache_hit.reset()
    proposal_cache_miss = Counter(
        'fragalysis_proposal_cache_miss',
        'Number of proposal cache misses',
    )
    proposal_cache_miss.reset()

    @staticmethod
    def new_tunnel():
        PrometheusMetrics.ssh_tunnels.inc()

    @staticmethod
    def failed_tunnel():
        PrometheusMetrics.ssh_tunnel_failures.inc()

    @staticmethod
    def new_ispyb_connection():
        PrometheusMetrics.ispyb_connections.inc()

    @staticmethod
    def new_ispyb_connection_attempt():
        PrometheusMetrics.ispyb_connection_attempts.inc()

    @staticmethod
    def failed_ispyb_connection():
        PrometheusMetrics.ispyb_connection_failures.inc()

    @staticmethod
    def new_proposal_cache_hit():
        PrometheusMetrics.proposal_cache_hit.inc()

    @staticmethod
    def new_proposal_cache_miss():
        PrometheusMetrics.proposal_cache_miss.inc()
