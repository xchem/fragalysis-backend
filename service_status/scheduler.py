import logging

from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.interval import IntervalTrigger

logger = logging.getLogger('service_status')

_state: dict = {}


def get_scheduler() -> BackgroundScheduler:
    if 'scheduler' not in _state:
        _state['scheduler'] = BackgroundScheduler(daemon=True)
    return _state['scheduler']


def add_service_job(service_name: str, func, interval_seconds: int) -> None:
    """Register (or replace) an interval job for a single service."""
    scheduler = get_scheduler()
    job_id = f"service_status.{service_name}"
    scheduler.add_job(
        func,
        trigger=IntervalTrigger(seconds=interval_seconds),
        id=job_id,
        name=service_name,
        replace_existing=True,
        max_instances=1,
        misfire_grace_time=10,
    )


def start() -> None:
    scheduler = get_scheduler()
    if not scheduler.running:
        logger.info('Starting APScheduler BackgroundScheduler')
        scheduler.start()


def shutdown(wait: bool = True) -> None:
    scheduler = _state.get('scheduler')
    if scheduler and scheduler.running:
        scheduler.shutdown(wait=wait)
