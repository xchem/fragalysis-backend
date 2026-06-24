"""LoggingAdapters used by various modules.
"""

import logging


class TaskLoggerAdapter(logging.LoggerAdapter):
    """
    An adapter used by modules normally used in Celery tasks. Given a dictionary that
    provides at least a 'task' (a UUID string). It also supports 'target', 'tas', and
    'username' keys.

    The caller can also provide a contextual 'marker', added to each logged entry to
    allow lines relating to a simply objective to be seen. For example,
    you might provide the marker `DOWNLOAD` and all the lines will then contain the
    word `DOWNLOAD`.

    The generated prefix string always contains the first 8 characters of the task
    ID (i.e. '54e3ea08-383a-4d96-9fc9-51948c8130b6' becomes '54e3ea08') and
    is printed as 'T.<task ID>'. The prefix also optionally contains
    't(<target>)', 'tas(<tas>)', and 'u(<username>)'. If the username exists but
    is blank, the username is replaced with '|anon|'.
    """

    def process(self, msg, kwargs):
        assert self.extra
        if 'task' in self.extra:
            rendered_extra: str = f"T.{str(self.extra['task'])[:8]}"
        else:
            rendered_extra = "T.000000"
        if 'marker' in self.extra:
            rendered_extra += f" {self.extra['marker']}"
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
