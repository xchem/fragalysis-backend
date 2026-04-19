"""LoggingAdapters used by various modules.
"""

import logging


class TaskLoggerAdapter(logging.LoggerAdapter):
    """
    An adapter used by modules normally used in Celery tasks. Given a dictionary that
    provides at least a 'task' (a UUID string). It also supports 'target', 'tas', and
    'username' keys.

    The generated string always contains '.T.<task ID>' that optionally also contains
    't(<target>)', 'tas(<tas>)', and 'u(<username>)'. If the username exists but
    is blank, the username is replaced with '|anon|'.
    """

    def process(self, msg, kwargs):
        assert self.extra
        rendered_extra: str = f"'.T.{self.extra['task']}"
        if 'target' in self.extra:
            rendered_extra += f" t({self.extra['target']})"
        if 'tas' in self.extra:
            rendered_extra += f" tas({self.extra['tas']})"
        if 'username' in self.extra:
            if self.extra['username']:
                rendered_extra += f" u({self.extra['username']})"
            else:
                rendered_extra += " u(|anon|)"
        return (rendered_extra, kwargs)
