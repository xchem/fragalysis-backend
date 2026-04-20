"""LoggingAdapters used by various modules.
"""

import logging


class TaskLoggerAdapter(logging.LoggerAdapter):
    """
    An adapter used by modules normally used in Celery tasks. Given a dictionary that
    provides at least a 'task' (a UUID string). It also supports 'target', 'tas', and
    'username' keys.

    The generated prefix string always contains the last 6 characters of the Task UUID
    string printed as '.T.<task ID>'. The string also optionally contains
    't(<target>)', 'tas(<tas>)', and 'u(<username>)'. If the username exists but
    is blank, the username is replaced with '|anon|'.
    """

    def process(self, msg, kwargs):
        assert self.extra
        if 'task' in self.extra and len(str(self.extra['task'])) > 6:
            rendered_extra: str = f"'.T.{str(self.extra['task'])[-6:]}"
        else:
            rendered_extra = ".T.000000"
        if 'target' in self.extra:
            rendered_extra += f" t({self.extra['target']})"
        if 'tas' in self.extra:
            rendered_extra += f" tas({self.extra['tas']})"
        if 'username' in self.extra:
            if self.extra['username']:
                rendered_extra += f" u({self.extra['username']})"
            else:
                rendered_extra += " u(|anon|)"
        return '%s %s' % (rendered_extra, msg), kwargs
