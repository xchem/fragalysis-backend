import logging

from django.db import models

from .managers import ServiceStateDataManager

# TODO: separate file
logger = logging.getLogger(__name__)


class ServiceState(models.Model):
    """Query result choices:
    - NOT_CONFIGURED
    - DEGRADED
    - OK
    - ERROR
    """

    state = models.TextField(null=False, primary_key=True)
    display_name = models.TextField(null=False)

    def is_success(self) -> bool:
        return self.state == 'OK'


class Service(models.Model):
    service = models.TextField(null=False, primary_key=True)
    display_name = models.TextField(null=False)
    frequency = models.PositiveSmallIntegerField(
        default=30, help_text='Ping frequency in seconds'
    )
    last_state = models.ForeignKey(
        ServiceState, null=False, on_delete=models.PROTECT, default='NOT_CONFIGURED'
    )
    last_states_of_same_type = models.IntegerField(null=False, default=0)
    last_query_time = models.DateTimeField(null=True)
    last_success = models.DateTimeField(null=True)
    last_failure = models.DateTimeField(null=True)

    objects = models.Manager()
    data_manager = ServiceStateDataManager()

    def __str__(self) -> str:
        return f"{self.service}"

    def __repr__(self) -> str:
        return "<Service %r %r>" % (self.service, self.last_result)
