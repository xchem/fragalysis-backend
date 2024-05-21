# import inspect
import logging

from django.apps import AppConfig
from django.conf import settings

logger = logging.getLogger(__name__)


class ServiceStatusConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'service_status'

    @staticmethod
    def start_service_queries():
        import service_status.services as services

        from .models import Service

        # from .models import ServiceState
        # from .services import State

        logger.warning('starting services')
        service_string = settings.ENABLE_SERVICE_STATUS
        requested_services = [k for k in service_string.split(":") if k != ""]

        # The next for block is supposed to check if the service check
        # test is written in services.py and if it's requested by env
        # variable, but something weird happens: The line
        # service = Service.objects.get(service=name)
        # triggers an error: relation Services doesn't exist. Weirdly,
        # the same call below, in next for block succeeds and even
        # more weirdly, if I move this block below the next one, it
        # still fails. As if the table exists and then gets deleted. I
        # have no idea what is happening and I cannot find any
        # meaningful difference between them.

        # test function exists but not requested, set not configured
        # for name, body in inspect.getmembers(services, inspect.isfunction):
        #     if (
        #         inspect.getsource(body).find('@shared_task')
        #         and inspect.getsource(body).find('@service_query')
        #         and name not in requested_services
        #     ):
        #         try:
        #             service = Service.objects.get(service=name)
        #             service.last_state = ServiceState.objects.get(
        #                 state=State.NOT_CONFIGURED
        #             )
        #             service.save()
        #         except Service.DoesNotExist:
        #             # populate if missing
        #             Service(
        #                 service=name,
        #                 display_name=name.capitalize().replace('_', ' '),
        #             ).save()

        # service is requested but query function doesn't exist, also add not configured
        for service in requested_services:
            logger.warning('configuring service: %s', service)
            task = getattr(services, service, None)
            if not task:
                logger.error(
                    'Task for service %s does not exist, test function not implemented',
                    service,
                )
                continue

            logger.warning('found task: %s', task)

            # by now it should be possible to launch a service query,
            # check if exists in db, add if not
            try:
                _ = Service.objects.get(service=service)
            except Service.DoesNotExist:
                # populate if missing
                # duplicate key error?? how??
                Service(
                    service=service,
                    display_name=service.capitalize().replace('_', ' '),
                ).save()

            task.delay(1)

    def ready(self):
        logger.warning('ready method running')
        self.start_service_queries()
