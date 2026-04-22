import concurrent.futures
import functools
import inspect
import logging
import os
from enum import Enum

from django.conf import settings
from django.utils import timezone

from .models import Service, ServiceState

# What's our HOSTNAME?
# If it's _SERVICE_CHECK_HOSTNAME then we start services, otherwise we don't
_HOSTNAME: str = os.environ.get('HOSTNAME', '')
_SERVICE_CHECK_HOSTNAME: str = 'stack-0'

logger = logging.getLogger('service_status')

# Soft timeout for individual service query calls
SERVICE_QUERY_TIMEOUT_S = 28

_executor = concurrent.futures.ThreadPoolExecutor(max_workers=10)


# this is a bit redundant because they're all in database, but it's
# convenient to have them here
class State(str, Enum):
    NOT_CONFIGURED = "NOT_CONFIGURED"
    DEGRADED = "DEGRADED"
    OK = "OK"
    ERROR = "ERROR"


def service_query(func):
    """Decorator function for service queries functions"""

    @functools.wraps(func)
    def wrapper_service_query(*args, **kwargs):  # pylint: disable=unused-argument
        service = Service.objects.get(service=func.__name__)

        # If the service has been disabled, skip the check
        if service.last_state_id == State.NOT_CONFIGURED:
            logger.debug('Service %s is NOT_CONFIGURED, skipping', func.__name__)
            return

        try:
            future = _executor.submit(func)
            state_pk = future.result(timeout=SERVICE_QUERY_TIMEOUT_S)
        except concurrent.futures.TimeoutError:
            logger.warning('Query time limit exceeded, setting result as DEGRADED')
            state_pk = State.DEGRADED

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

        service.total_queries = service.total_queries + 1

        # unexplored possibility to adjust ping times if necessary

        service.save()

    wrapper_service_query.is_service_query = True
    return wrapper_service_query


def init_services():
    logger.debug('+ init_services')

    # Do nothing if we're not the service check Pod.
    # Only one Pod needs to check the service status.
    if _HOSTNAME != _SERVICE_CHECK_HOSTNAME:
        logger.warning(
            'This host (%s) is not the service check host (%s) - skipping initialisation',
            _HOSTNAME,
            _SERVICE_CHECK_HOSTNAME,
        )
        return

    service_string = settings.ENABLE_SERVICE_STATUS
    requested_services = [k for k in service_string.split(":") if k != ""]

    import service_status.services as services_module

    # gather all service functions from services.py and make sure they're
    # in db; identified by the _is_service_query marker set by @service_query
    defined = []
    for name, obj in inspect.getmembers(services_module, inspect.isfunction):
        if getattr(obj, 'is_service_query', False):
            defined.append(name)
            # ensure all defined services are in db
            try:
                service = Service.objects.get(service=name)
            except Service.DoesNotExist:
                # add if missing
                docs = inspect.getdoc(obj)
                display_name = docs.splitlines()[0] if docs else ''
                Service(
                    service=name,
                    # use first line of docstring as user-friendly name
                    display_name=display_name,
                ).save()

    # clear up db of those that are not defined
    for service in Service.objects.all():
        if service.service not in defined:
            service.delete()

    # mark those not requested as NOT_CONFIGURED
    for service in Service.objects.all():
        if service.service not in requested_services:
            service.last_state = ServiceState.objects.get(state=State.NOT_CONFIGURED)
            service.save()

    # Now start the scheduler and add jobs for requested services
    if settings.SERVICE_STATUS_SCHEDULER_ENABLED:
        import service_status.scheduler as scheduler_module

        scheduler_module.start()
        for s in requested_services:
            logger.debug('Trying to schedule service: %s', s)
            try:
                service = Service.objects.get(service=s)
            except Service.DoesNotExist:
                logger.error(
                    'Service %s requested but function missing in services.py',
                    s,
                )
                continue

            func = getattr(services_module, service.service)
            logger.debug(
                'Adding scheduler job for: %s (every %ds)', s, service.frequency
            )
            scheduler_module.add_service_job(service.service, func, service.frequency)


def services(enable=(), disable=()):
    logger.debug('+ services')

    if enable is None:
        enable = []
    if disable is None:
        disable = []

    to_enable = set(enable).difference(set(disable))
    to_disable = set(disable).difference(set(enable))
    confusables = set(disable).intersection(set(enable))

    # at this point, all the services must be started and in db
    for name in to_enable:
        try:
            service = Service.objects.get(service=name)
        except Service.DoesNotExist:
            logger.error('Unknown service: %s', name)
            continue

        if service.last_state_id == State.NOT_CONFIGURED:
            service.last_state = ServiceState.objects.get(state=State.DEGRADED)
            service.save()
            logger.info('Enabled service %s', name)
        else:
            logger.info('Service %s is already enabled', name)

    for name in to_disable:
        try:
            service = Service.objects.get(service=name)
        except Service.DoesNotExist:
            logger.error('Unknown service: %s', name)
            continue

        service.last_state = ServiceState.objects.get(state=State.NOT_CONFIGURED)
        service.save()
        logger.info('Disabled service %s', name)

    # service name in both enable and disable: toggle based on current state
    for name in confusables:
        try:
            service = Service.objects.get(service=name)
        except Service.DoesNotExist:
            logger.error('Unknown service: %s', name)
            continue

        if service.last_state_id == State.NOT_CONFIGURED:
            service.last_state = ServiceState.objects.get(state=State.DEGRADED)
            service.save()
            logger.info('Enabled service %s (was NOT_CONFIGURED)', name)
        else:
            service.last_state = ServiceState.objects.get(state=State.NOT_CONFIGURED)
            service.save()
            logger.info('Disabled service %s', name)
