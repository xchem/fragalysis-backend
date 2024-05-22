import functools
import inspect
import logging
from enum import Enum

from celery.exceptions import SoftTimeLimitExceeded
from django.conf import settings
from django.utils import timezone

from .models import Service, ServiceState

logger = logging.getLogger('service_status')


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
        import service_status.services as services_module

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

        # get the task function from this module
        task = getattr(services_module, func.__name__)
        task.delay(service.frequency)

    return wrapper_service_query


def init_services():
    logger.debug('+ init_services')
    service_string = settings.ENABLE_SERVICE_STATUS
    requested_services = [k for k in service_string.split(":") if k != ""]

    import service_status.services as services_module

    # gather all test functions from services.py and make sure they're
    # in db
    defined = []
    for name, body in inspect.getmembers(services_module):
        # doesn't seem to have a neat way to test if object is task,
        # have to check them manually
        try:
            src = inspect.getsource(body)
        except TypeError:
            # uninteresting propery
            continue

        if src.find('@shared_task') >= 0 and src.find('@service_query') >= 0:
            defined.append(name)
            # ensure all defined services are in db
            try:
                service = Service.objects.get(service=name)
            except Service.DoesNotExist:
                # add if missing
                docs = inspect.getdoc(body)
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

    # and launch the rest
    # TODO: this could potentially be an actual check if beat is running
    if not settings.CELERY_TASK_ALWAYS_EAGER:
        for s in requested_services:
            logger.debug('trying to launch service: %s', s)
            try:
                service = Service.objects.get(service=s)
            except Service.DoesNotExist:
                logger.error(
                    'Service %s requested but test function missing in services.py',
                    s,
                )
                continue

            # launch query task
            task = getattr(services_module, service.service)
            logger.debug('trying to launch task: %s', task)
            task.delay(1)


def services(enable=(), disable=()):
    logger.debug('+ init_services')
    import service_status.services as services_module

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

        task = getattr(services_module, service.service)
        task.delay(1)
        logger.info('Starting service query %s', name)

    # TODO: this bit isn't working really
    for name in to_disable:
        try:
            service = Service.objects.get(service=name)
        except Service.DoesNotExist:
            logger.error('Unknown service: %s', name)
            continue

        task = getattr(services_module, service.service)

    for name in confusables:
        try:
            service = Service.objects.get(service=name)
        except Service.DoesNotExist:
            logger.error('Unknown service: %s', name)
            continue

        task = getattr(services_module, service.service)
