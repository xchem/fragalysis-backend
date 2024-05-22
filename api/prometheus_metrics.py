"""Prometheus metrics used by the fragalysis API module.
"""
from django.conf import settings
from prometheus_client import Counter

_NAMESPACE: str = settings.OUR_KUBERNETES_NAMESPACE or 'unknown'


class PrometheusMetrics:
    ssh_tunnels = Counter(
        'fragalysis_ssh_tunnels',
        'Total number of SSH tunnels created',
        labelnames=['k8s_namespace'],
    )
    ssh_tunnel_failures = Counter(
        'fragalysis_ssh_tunnel_failures',
        'Total number of SSH tunnel failures',
        labelnames=['k8s_namespace'],
    )
    ispyb_connections = Counter(
        'fragalysis_ispyb_connections',
        'Total number of ISpyB connections',
        labelnames=['k8s_namespace'],
    )
    ispyb_connection_failures = Counter(
        'fragalysis_ispyb_connection_failures',
        'Total number of ISpyB connection failures',
        labelnames=['k8s_namespace'],
    )

    @staticmethod
    def new_tunnel():
        PrometheusMetrics.ssh_tunnels.labels(_NAMESPACE).inc()

    @staticmethod
    def failed_tunnel():
        PrometheusMetrics.ssh_tunnel_failures.labels(_NAMESPACE).inc()

    @staticmethod
    def new_ispyb_connection():
        PrometheusMetrics.ispyb_connections.labels(_NAMESPACE).inc()

    @staticmethod
    def failed_ispyb_connection():
        PrometheusMetrics.ispyb_connection_failures.labels(_NAMESPACE).inc()
