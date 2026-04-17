"""LoggingAdapters used by various modules.
"""

import logging


class TaskLoggerAdapter(logging.LoggerAdapter):
    """
    An adapter used by modules normally used in Celery tasks,
    which expects a dictionary that provide 'task', 'target' and 'tas'.
    """

    def process(self, msg, kwargs):
        assert self.extra
        return (
            '/%s/ TAS(%s) Tgt(%s) %s'
            % (self.extra['task'], self.extra['tas'], self.extra['target'], msg),
            kwargs,
        )
