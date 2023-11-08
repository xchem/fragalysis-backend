import os
import logging
from enum import Enum

import asyncio
import functools
import requests

from pydiscourse import DiscourseClient

from viewer.squonk2_agent import get_squonk2_agent
from frag.utils.network_utils import get_driver
from api import security

logger = logging.getLogger(__name__)

DELAY = 10


class State(str, Enum):
    NOT_CONFIGURED = "NOT_CONFIGURED"
    DEGRADED = "DEGRADED"
    OK = "OK"
    TIMEOUT = "TIMEOUT"
    ERROR = "ERROR"


class Service(str, Enum):
    ISPYB = "ispyb"
    DISCOURSE = "discourse"
    SQUONK = "squonk"
    FRAG = "fragmentation_graph"
    KEYCLOAK = "keycloak"


# called from the outside
def get_service_state(services):
    return asyncio.run(service_queries(services))


async def service_queries(services):
    """Chain the requested service queries"""
    logger.debug("service query called")
    coroutines = []
    if Service.ISPYB in services:
        coroutines.append(
            ispyb(
                Service.ISPYB,
                "Access control (ISPyB)",
                ispyb_host="ISPYB_HOST",
            )
        )

    if Service.SQUONK in services:
        coroutines.append(
            squonk(
                Service.SQUONK,
                "Squonk",
                squonk_pwd="SQUONK2_ORG_OWNER_PASSWORD",
            )
        )

    if Service.FRAG in services:
        coroutines.append(
            fragmentation_graph(
                Service.FRAG,
                "Fragmentation graph",
                url="NEO4J_BOLT_URL",
            )
        )

    if Service.DISCOURSE in services:
        coroutines.append(
            discourse(
                Service.DISCOURSE,
                "Discourse",
                key="DISCOURSE_API_KEY",
                url="DISCOURSE_HOST",
                user="DISCOURSE_USER",
            )
        )

    if Service.KEYCLOAK in services:
        coroutines.append(
            keycloak(
                Service.KEYCLOAK,
                "Keycloak",
                url="OIDC_KEYCLOAK_REALM",
                secret="OIDC_RP_CLIENT_SECRET",
            )
        )

    result = await asyncio.gather(*coroutines)
    logger.debug("service query result: %s", result)
    return result


def service_query(func):
    """Decorator function for service queries functions"""

    @functools.wraps(func)
    async def wrapper_service_query(*args, **kwargs):
        logger.debug("+ wrapper_service_query")
        logger.debug("args passed: %s", args)
        logger.debug("kwargs passed: %s", kwargs)
        logger.debug("function: %s", func.__name__)

        state = State.NOT_CONFIGURED
        envs = [os.environ.get(k, None) for k in kwargs.values()]
        # env variables come in kwargs, all must be defined
        if all(envs):
            state = State.DEGRADED
            loop = asyncio.get_running_loop()
            try:
                async with asyncio.timeout(DELAY):
                    future = loop.run_in_executor(
                        None, functools.partial(func, *args, **kwargs)
                    )
                    logger.debug("future: %s", future)
                    result = await future
                    logger.debug("result: %s", result)
                    if result:
                        state = State.OK

            except TimeoutError:
                logger.error("Service query '%s' timed out", func.__name__)
                state = State.TIMEOUT
            except Exception as exc:
                # unknown error with the query
                logger.exception(exc, exc_info=True)
                state = State.ERROR

        # name and ID are 1nd and 0th params respectively.
        # alternative solution for this would be to return just a
        # state and have the service_queries() map the results to the
        # correct values
        return {"id": args[0], "name": args[1], "state": state}

    return wrapper_service_query


@service_query
async def ispyb(func_id, name, ispyb_host=None):
    logger.debug("+ ispyb")
    return security.get_conn()


@service_query
def discourse(func_id, name, key=None, url=None, user=None):
    logger.debug("+ discourse")
    client = DiscourseClient(
        os.environ.get(url, None),
        api_username=os.environ.get(user, None),
        api_key=os.environ.get(key, None),
    )
    # TODO: some action on client?
    return client


@service_query
def squonk(func_id, name, squonk_pwd=None):
    logger.debug("+ squonk")
    return get_squonk2_agent().configured().success


@service_query
def fragmentation_graph(func_id, name, url=None):
    logger.debug("+ fragmentation_graph")
    # graph_driver = get_driver(url='neo4j', neo4j_auth='neo4j/password')
    graph_driver = get_driver()
    with graph_driver.session() as session:
        try:
            _ = session.run("match (n) return count (n);")
            return True
        except ValueError:
            # service isn't running
            return False


@service_query
def keycloak(func_id, name, url=None, secret=None):
    logger.debug("+ keycloak")
    keycloak_realm = os.environ.get(url, None)
    response = requests.get(keycloak_realm)
    logger.debug("keycloak response: %s", response)
    return response.ok
