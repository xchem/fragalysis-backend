import logging
import time
from random import random

import requests
from celery import shared_task
from django.conf import settings
from frag.utils.network_utils import get_driver
from pydiscourse import DiscourseClient

from api.security import ping_configured_connector
from viewer.squonk2_agent import get_squonk2_agent

from .utils import State, service_query

logger = logging.getLogger('service_status')


# Default timeout for any request calls
# Used for keycloak atm.
REQUEST_TIMEOUT_S = 5

# Service query timeout
SERVICE_QUERY_TIMEOUT_S = 28


# service status test functions
# NB! first line of docstring is used as a display name


# @shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
# @service_query
def test_query() -> str:
    """A dumb little test query.

    For testing.
    """
    logger.debug('+ test_query')
    state = State.DEGRADED
    time.sleep(3)
    if random() > 0.2:
        state = State.ERROR
        if random() > 0.2:
            state = State.OK
        else:
            state = State.ERROR

    logger.debug('end state: %s', state)
    return state.name


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def ispyb() -> str:
    """Access control (ISPyB)"""
    logger.debug("+ ispyb")
    return State.OK if ping_configured_connector() else State.DEGRADED


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def discourse() -> str:
    """Discourse"""
    logger.debug("+ discourse")
    # Discourse is "unconfigured" if there is no API key
    if not any(
        [settings.DISCOURSE_API_KEY, settings.DISCOURSE_HOST, settings.DISCOURSE_USER]
    ):
        return State.NOT_CONFIGURED
    client = DiscourseClient(
        settings.DISCOURSE_HOST,
        api_username=settings.DISCOURSE_USER,
        api_key=settings.DISCOURSE_API_KEY,
    )
    # TODO: some action on client?
    return State.DEGRADED if client is None else State.OK


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def squonk() -> str:
    """Squonk"""
    logger.debug("+ squonk")
    return State.OK if get_squonk2_agent().configured().success else State.DEGRADED


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def fragmentation_graph() -> str:
    """Fragmentation graph"""
    logger.debug("+ fragmentation_graph")
    graph_driver = get_driver(url=settings.NEO4J_QUERY, neo4j_auth=settings.NEO4J_AUTH)
    with graph_driver.session() as session:
        try:
            _ = session.run("match (n) return count (n);")
            return State.OK
        except ValueError:
            # service isn't running
            return State.DEGRADED


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def keycloak() -> str:
    """Keycloak"""
    logger.debug("+ keycloak")
    # Keycloak is "unconfigured" if there is no realm URL
    keycloak_realm = settings.OIDC_KEYCLOAK_REALM
    if not keycloak_realm:
        return State.NOT_CONFIGURED
    response = requests.get(keycloak_realm, timeout=REQUEST_TIMEOUT_S)
    logger.debug("keycloak response: %s", response)
    return State.OK if response.ok else State.DEGRADED
