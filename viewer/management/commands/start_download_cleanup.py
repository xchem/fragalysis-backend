"""Using the service status scheduler we install tasks that run regularly to
clean-up "out od date" DownloadLinks records, and old tmp/ media files.

A soft removal marks records as "expired" and a hard removal removes the underlying
files for "expired" records that have been expired for a significant time. The
"lost task" removal expires in-progress records whose Celery task has crashed or
been lost (so they will never finish). The tmp/ expiry deletes media storage
objects in tmp/ that are older than MEDIA_TMP_MAX_AGE_H hours.
"""
import signal

from django.conf import settings
from django.core.management.base import BaseCommand

import service_status.scheduler as scheduler_module
from viewer.download_structures import (
    expire_lost_download_records,
    hard_erase_out_of_date_download_records,
    soft_erase_out_of_date_download_records,
)
from viewer.media_cleanup import expire_tmp_files

_INTERVAL_S: int = settings.DOWNLOAD_CLEANUP_INTERVAL_M * 60


class Command(BaseCommand):
    help = "Periodically erase out-of-date download records and tmp/ media files (runs until terminated)"

    def handle(self, *args, **kwargs):
        del args, kwargs

        def _shutdown(signum, frame):  # pylint: disable=unused-argument
            self.stdout.write('Shutting down download cleanup scheduler...')
            scheduler_module.shutdown()

        signal.signal(signal.SIGTERM, _shutdown)
        signal.signal(signal.SIGINT, _shutdown)

        scheduler_module.add_service_job(
            'download_cleanup.soft_erase',
            soft_erase_out_of_date_download_records,
            _INTERVAL_S,
        )
        scheduler_module.add_service_job(
            'download_cleanup.hard_erase',
            hard_erase_out_of_date_download_records,
            _INTERVAL_S,
        )
        scheduler_module.add_service_job(
            'download_cleanup.expire_lost',
            expire_lost_download_records,
            _INTERVAL_S,
        )
        scheduler_module.add_service_job(
            'media_cleanup.expire_tmp',
            expire_tmp_files,
            _INTERVAL_S,
        )
        scheduler_module.start()

        self.stdout.write('Download cleanup scheduler started.')
        signal.pause()
