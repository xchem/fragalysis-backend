import os

_CONCURRENCY: int = int(os.environ.get("STACK_CONCURRENCY", "4"))

logconfig = "/code/gunicorn.log.conf"
bind = "unix:django_app.sock"
daemon = True
timeout = 3_000
workers = _CONCURRENCY
disable_redirect_access_to_syslog = True
