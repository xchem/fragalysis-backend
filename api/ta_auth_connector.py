"""A module that provides simplified request access to the TA Authenticator service.
Provides the following functions, that access the authenticator Pod: -

- get_auth_version()
- get_auth_ping()
- get_auth_target_access(username)
"""

import logging
from dataclasses import dataclass
from urllib.parse import quote

import requests
from django.conf import settings

_URL_TIMEOUT: int = 3
_URL_BASE: str = f'http://{settings.TA_AUTH_SERVICE}'
_QUERY_HEADERS: dict[str, str] = {'X-TAAQueryKey': settings.TA_AUTH_QUERY_KEY}

logger: logging.Logger = logging.getLogger(__name__)


@dataclass
class TasAuthVersionGetResponse:
    """The TA authenticator version response (including the base URL)."""

    version: str
    kind: str
    name: str
    location: str = _URL_BASE


@dataclass
class TasAuthPingGetResponse:
    """The TA authenticator ping response."""

    ping: str


def get_auth_version() -> TasAuthVersionGetResponse:
    """Returns the version reported by the TA authentication service."""
    if not settings.TA_AUTH_SERVICE:
        return TasAuthVersionGetResponse(
            version='', kind='SERVICE_NOT_PRESENT', name=''
        )
    if not settings.TA_AUTH_QUERY_KEY:
        return TasAuthVersionGetResponse(
            version='', kind='SERVICE_QUERY_KEY_NOT_DEFINED', name=''
        )

    url: str = f'{_URL_BASE}/version/'
    resp: requests.Response | None = None
    try:
        resp = requests.get(url, timeout=_URL_TIMEOUT)
    except requests.exceptions.RequestException as r_ex:  # pylint: disable=broad-except
        logger.error('TA:GET:%s RequestException (%s)', url, r_ex)
    except Exception as ex:  # pylint: disable=broad-exception-caught
        logger.error('TA:GET:%s Exception (%s)', url, ex)

    if resp is None:
        logger.warning('TA:GET:%s (no response)', url)
        return TasAuthVersionGetResponse(
            version='', kind='ERROR_INTERNAL', name='Null response'
        )
    elif resp.status_code not in (200,):
        logger.warning('TA:GET:%s [%s] (status not 200)', url, resp.status_code)
        return TasAuthVersionGetResponse(
            version='', kind='ERROR_INTERNAL', name='(status not 200)'
        )
    elif 'application/json' not in resp.headers.get('Content-Type', ''):
        logger.warning('TA:GET:%s (empty response)', url)
        return TasAuthVersionGetResponse(
            version='', kind='ERROR_INTERNAL', name='(empty response)'
        )
    elif 'version' not in resp.json():
        logger.warning('TA:GET:%s (no version property)', url)
        return TasAuthVersionGetResponse(
            version='', kind='ERROR_INTERNAL', name='(no version property)'
        )

    logger.info('TA:GET:%s [%s]', url, resp.json())

    return TasAuthVersionGetResponse(
        version=resp.json()['version'],
        kind=resp.json()['kind'],
        name=resp.json()['name'],
    )


def get_auth_ping() -> TasAuthPingGetResponse:
    """Returns the ping reported by the TA authentication service."""
    if not settings.TA_AUTH_SERVICE:
        return TasAuthPingGetResponse(ping='SERVICE_NOT_PRESENT')

    url: str = f'{_URL_BASE}/ping/'
    resp: requests.Response | None = None
    try:
        resp = requests.get(url, timeout=_URL_TIMEOUT)
    except requests.exceptions.RequestException as r_ex:  # pylint: disable=broad-except
        logger.error('TA:GET:%s RequestException (%s)', url, r_ex)
    except Exception as ex:  # pylint: disable=broad-exception-caught
        logger.error('TA:GET:%s Exception (%s)', url, ex)

    if resp is None:
        logger.warning('TA:GET:%s (no response)', url)
        return TasAuthPingGetResponse('PING response was null')
    elif resp.status_code not in (200,):
        logger.warning('TA:GET:%s [%s] (status not 200)', url, resp.status_code)
        return TasAuthPingGetResponse('PING response status not 200')
    elif 'application/json' not in resp.headers.get('Content-Type', ''):
        logger.warning('TA:GET:%s (empty response)', url)
        return TasAuthPingGetResponse('PING response was empty')
    elif 'ping' not in resp.json():
        logger.warning('TA:GET:%s (no ping property)', url)
        return TasAuthPingGetResponse('PING response has no ping property')

    ping: str = resp.json()['ping']
    logger.info('TA:GET:%s [%s]', url, ping)

    return TasAuthPingGetResponse(ping=ping)


def get_auth_target_access(username: str) -> set[str]:
    """Returns the set of target access strings a user is entitled to
    as reported by the TA authentication service."""
    assert username

    empty_target_access: set[str] = set()

    if not settings.TA_AUTH_QUERY_KEY:
        logger.debug('Skipping query - query key is not set (TA_AUTH_QUERY_KEY)')
        return empty_target_access

    url: str = f'{_URL_BASE}/target-access/{quote(username)}'
    resp: requests.Response | None = None
    try:
        resp = requests.get(url, headers=_QUERY_HEADERS, timeout=_URL_TIMEOUT)
    except requests.exceptions.RequestException as ex:  # pylint: disable=broad-except
        logger.error('TA:GET:%s RequestException (%s)', url, ex)
    except Exception as ex:  # pylint: disable=broad-exception-caught
        logger.error('TA:GET:%s Exception (%s)', url, ex)

    if resp is None:
        logger.warning('TA:GET:%s (no response)', url)
        return empty_target_access
    if resp.status_code not in (200,):
        logger.warning('TA:GET:%s [%s] (status not 200)', url, resp.status_code)
        return empty_target_access
    elif 'application/json' not in resp.headers.get('Content-Type', ''):
        logger.warning('TA:GET:%s (empty response)', url)
        return empty_target_access
    elif 'count' not in resp.json():
        logger.warning('TA:GET:%s (no count)', url)
        return empty_target_access
    elif 'target_access' not in resp.json():
        logger.warning('TA:GET:%s (no target_access)', url)
        return empty_target_access

    target_access: set[str] = set(resp.json()['target_access'])
    logger.info('TA:GET:%s (got %d for %s)', url, len(target_access), username)

    return target_access
