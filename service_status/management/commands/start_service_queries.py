import os
import signal

from django.core.management.base import BaseCommand

import service_status.scheduler as scheduler_module
from service_status.utils import init_services

_HOSTNAME: str = os.environ.get('HOSTNAME', '')
_SERVICE_CHECK_HOSTNAME: str = 'stack-0'


class Command(BaseCommand):
    help = "Activate service health check queries (runs until terminated)"

    def handle(self, *args, **kwargs):
        # Unused args
        del args, kwargs

        if _HOSTNAME != _SERVICE_CHECK_HOSTNAME:
            self.stdout.write(
                f'This host ({_HOSTNAME}) is not the service check host'
                f' ({_SERVICE_CHECK_HOSTNAME}) - doing nothing'
            )
            return

        def _shutdown(signum, frame):  # pylint: disable=unused-argument
            self.stdout.write('Shutting down service scheduler...')
            scheduler_module.shutdown()

        signal.signal(signal.SIGTERM, _shutdown)
        signal.signal(signal.SIGINT, _shutdown)

        init_services()
        self.stdout.write('Service scheduler started.')
        signal.pause()
