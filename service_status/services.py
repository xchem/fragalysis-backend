import functools
import logging
import os
import sys
import time
from enum import Enum
from random import random

import requests
from celery import shared_task
from celery.exceptions import SoftTimeLimitExceeded
from django.conf import settings
from django.utils import timezone
from frag.utils.network_utils import get_driver
from pydiscourse import DiscourseClient

from api.security import ping_configured_connector
from viewer.squonk2_agent import get_squonk2_agent

from .models import Service, ServiceState

logger = logging.getLogger('service_status')


# Default timeout for any request calls
# Used for keycloak atm.
REQUEST_TIMEOUT_S = 5

# Service query timeout
SERVICE_QUERY_TIMEOUT_S = 10


class State(str, Enum):
    NOT_CONFIGURED = "NOT_CONFIGURED"
    DEGRADED = "DEGRADED"
    OK = "OK"
    ERROR = "ERROR"


def service_query(func):
    """Decorator function for service queries functions"""

    @functools.wraps(func)
    def wrapper_service_query(*args, **kwargs):
        logger.debug("+ wrapper_service_query")
        logger.debug("args passed: %s", args)
        logger.debug("kwargs passed: %s", kwargs)
        logger.debug("function: %s", func.__name__)

        try:
            state_pk = func()
        except SoftTimeLimitExceeded:
            logger.warning('Query time limit exceeded, setting result as DEGRADED')
            state_pk = State.DEGRADED

        service = Service.objects.get(service=func.__name__)

        state = ServiceState.objects.get(state=state_pk)
        if service.last_state == state:
            service.last_states_of_same_type = service.last_states_of_same_type + 1
        else:
            service.last_states_of_same_type = 0

        service.last_state = state
        timestamp = timezone.now()
        service.last_query_time = timestamp
        if state.is_success():
            service.last_success = timestamp
        else:
            service.last_failure = timestamp

        # unexplored possibility to adjust ping times if necessary

        service.save()

        task = getattr(sys.modules[__name__], func.__name__)
        task.delay(service.frequency)

    return wrapper_service_query


# service status test functions


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def test_query():
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
def ispyb() -> bool:
    logger.debug("+ ispyb")
    return ping_configured_connector()


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def discourse() -> bool:
    logger.debug("+ discourse")
    # Discourse is "unconfigured" if there is no API key
    if not settings.DISCOURSE_API_KEY:
        return False
    client = DiscourseClient(
        os.environ.get(settings.DISCOURSE_HOST, None),
        api_username=os.environ.get(settings.DISCOURSE_USER, None),
        api_key=os.environ.get(settings.DISCOURSE_API_KEY, None),
    )
    # TODO: some action on client?
    return client != None


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def squonk() -> bool:
    logger.debug("+ squonk")
    return get_squonk2_agent().configured().success


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def fragmentation_graph() -> bool:
    logger.debug("+ fragmentation_graph")
    graph_driver = get_driver(
        url=settings.NEO4J_LOCATION, neo4j_auth=settings.NEO4J_AUTH
    )
    with graph_driver.session() as session:
        try:
            _ = session.run("match (n) return count (n);")
            return True
        except ValueError:
            # service isn't running
            return False


@shared_task(soft_time_limit=SERVICE_QUERY_TIMEOUT_S)
@service_query
def keycloak() -> bool:
    logger.debug("+ keycloak")
    # Keycloak is "unconfigured" if there is no realm URL
    keycloak_realm = os.environ.get(settings.OIDC_KEYCLOAK_REALM, None)
    if not keycloak_realm:
        return False
    response = requests.get(keycloak_realm, timeout=REQUEST_TIMEOUT_S)
    logger.debug("keycloak response: %s", response)
    return response.ok
