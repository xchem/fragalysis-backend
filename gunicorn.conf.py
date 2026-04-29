import logging
import os
import sys

_CONCURRENCY: int = int(os.environ.get("STACK_CONCURRENCY", "4"))

_GUNICORN_LOGGING_DIR: str = os.environ.get("GUNICORN_LOGGING_DIR", "/code/logs")
_GUNICORN_ROOT_LOGGING_LEVEL: str = os.environ.get(
    "GUNICORN_ROOT_LOGGING_LEVEL", "INFO"
).upper()
_GUNICORN_ERROR_LOGGING_LEVEL: str = os.environ.get(
    "GUNICORN_ERROR_LOGGING_LEVEL", "INFO"
).upper()
_GUNICORN_ACCESS_LOGGING_LEVEL: str = os.environ.get(
    "GUNICORN_ACCESS_LOGGING_LEVEL", "INFO"
).upper()
_GUNICORN_LOGGING_BACKUP_COUNT: int = int(
    os.environ.get("GUNICORN_LOGGING_BACKUP_COUNT", "7")
)

bind = "unix:django_app.sock"
daemon = True
timeout = 3_000
workers = _CONCURRENCY
disable_redirect_access_to_syslog = True


class LoggingPrometheusFilter(logging.Filter):
    def filter(self, record):
        return "GET /metrics" not in record.msg


logconfig_dict = {
    "version": 1,
    "disable_existing_loggers": True,
    "filters": {
        "prometheus_filter": {
            "()": LoggingPrometheusFilter,
        },
    },
    "formatters": {
        "generic": {
            "format": "%(asctime)s [%(process)d] [%(levelname)s] # %(message)s",
            "datefmt": "%Y-%m-%dT%H:%M:%S%z",
        },
        "access": {
            "format": "%(message)s",
        },
    },
    "handlers": {
        "console": {
            "level": _GUNICORN_ROOT_LOGGING_LEVEL,
            "class": "logging.StreamHandler",
            "stream": sys.stdout,
            "formatter": "generic",
        },
        "error_file": {
            "level": _GUNICORN_ERROR_LOGGING_LEVEL,
            "class": "logging.handlers.TimedRotatingFileHandler",
            "formatter": "generic",
            "filename": f"{_GUNICORN_LOGGING_DIR}/gunicorn-error.log",
            "when": "midnight",
            "interval": 1,
            "backupCount": _GUNICORN_LOGGING_BACKUP_COUNT,
            "encoding": "utf8",
        },
        "access_file": {
            "level": _GUNICORN_ACCESS_LOGGING_LEVEL,
            "class": "logging.handlers.TimedRotatingFileHandler",
            "formatter": "access",
            "filename": f"{_GUNICORN_LOGGING_DIR}/gunicorn-access.log",
            "when": "midnight",
            "interval": 1,
            "backupCount": _GUNICORN_LOGGING_BACKUP_COUNT,
            "encoding": "utf8",
            "filters": ["prometheus_filter"],
        },
    },
    "root": {
        "level": _GUNICORN_ROOT_LOGGING_LEVEL,
        "handlers": ["error_file"],
    },
    "loggers": {
        "root": {
            "level": _GUNICORN_ROOT_LOGGING_LEVEL,
            "handlers": ["error_file"],
        },
        "gunicorn.error": {
            "level": _GUNICORN_ERROR_LOGGING_LEVEL,
            "handlers": ["error_file"],
            "propagate": True,
        },
        "gunicorn.access": {
            "level": _GUNICORN_ACCESS_LOGGING_LEVEL,
            "handlers": ["access_file"],
            "propagate": False,
        },
    },
}
