"""LoggingAdapters used by various modules.
"""

import logging


class TaskLoggerAdapter(logging.LoggerAdapter):
    """
    An adapter used by modules normally used in Celery tasks. Given a dictionary that
    provides at least a 'task' (a UUID string). It also supports 'target', 'tas', and
    'username' keys.

    The generated string always contains '/<task ID>/' that optionally also contains
    'T(<target>)', 'TAS(<tas>)', and 'U(<username>)'.
    """

    def process(self, msg, kwargs):
        assert self.extra
        rendered_extra: str = f"'/{self.extra['task']}/"
        if 'target' in self.extra:
            rendered_extra += f" T({self.extra['target']})"
        if 'tas' in self.extra:
            rendered_extra += f" TAS({self.extra['tas']})"
        if 'username' in self.extra:
            rendered_extra += f" U({self.extra['username']})"
        return (rendered_extra, kwargs)
