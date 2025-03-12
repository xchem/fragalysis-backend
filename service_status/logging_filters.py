import logging

from django.conf import settings

_SERVICE_MODULE = "service_status"
_SERVICE_LOG_LEVEL = getattr(logging, settings.SERVICE_STATUS_LOGLEVEL, logging.WARNING)


# Ok, this took a while to figure out, so here are some notes for
# future me. This filter intercepts log messages and suppresses those
# associated with service_status module. They have two origins:
# 1) the module itself, the probe functions and all the support system
# 2) celery messages regarding task operations

# type 1 messages are intercepted with record.name attribute, this
# corresponds to service module name. In celery messages, the name
# attribute can have multiple values, so it's looking for the log
# message metadata in record.args, the name attribute there stores the
# function name that's generateing the message.

# Explored couple of avenues that turned out to be dead
# end/unnecessary for now: using the celery.current_task to get the
# task information. I added a flag to the functions in the
# @service_query decorator and checked the value here. This was able
# to identify the module functions but missed the celery messages
# which was half the problem. Another option would be to inspect the
# function body directly and check for the presence of decorators
# (like in service init function). Like said, this seems to be
# unnecessary now, but if more granularity is needed, this may be
# worth resucitating.


class SuppressServiceQueryTasksFilter(logging.Filter):
    def filter(self, record):
        if record.name == _SERVICE_MODULE:
            return record.levelno >= _SERVICE_LOG_LEVEL
        else:
            if isinstance(record.args, dict):
                name = record.args.get("name", "")
                if name.startswith(_SERVICE_MODULE):
                    return record.levelno >= _SERVICE_LOG_LEVEL

        return True  # let other messages through normally
