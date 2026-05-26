# pylint: disable=wildcard-import,unused-wildcard-import
"""Django settings used when running the test suite under pytest.

This module imports everything from the normal application settings and then
applies test-only overrides. It is selected via ``DJANGO_SETTINGS_MODULE`` in
``[tool.pytest.ini_options]`` (see ``pyproject.toml``).

The database connection details still come from the environment
(``POSTGRESQL_*``), which ``docker-compose.test.yml`` provides, because the
schema requires a real PostgreSQL instance (pgvector, simple_history, etc.) -
SQLite cannot stand in for it.
"""
from fragalysis.settings import *  # noqa: F401,F403

# Run as a non-production deployment. Production mode is intentionally
# stricter (see api.utils.deployment_mode_is_production); tests want the
# relaxed behaviour so they can exercise public/open-proposal code paths.
DEPLOYMENT_MODE = "DEVELOPMENT"

# Execute Celery tasks synchronously, in-process, and surface their
# exceptions to the caller rather than swallowing them.
CELERY_TASK_ALWAYS_EAGER = True
CELERY_TASK_EAGER_PROPAGATES = True

# By default the security layer resolves a user's proposals from Django
# (Project.user_id membership). Tests that want the external TA-authenticator
# path set settings.TA_AUTH_SERVICE themselves (e.g. via the pytest-django
# `settings` fixture) and patch api.ta_auth_connector.get_auth_target_access.
TA_AUTH_SERVICE = ""
TA_AUTH_QUERY_KEY = ""

# Use the plain PostgreSQL backend for tests instead of the django_prometheus
# wrapper. The wrapper registers process-wide Prometheus collectors at import
# time, which clashes ("Duplicated timeseries in CollectorRegistry") when
# settings are imported repeatedly under pytest. The schema is identical.
DATABASES["default"]["ENGINE"] = "django.db.backends.postgresql"  # noqa: F405

# Fast, insecure password hashing - we never need real hashing strength in
# tests and this noticeably speeds up user creation.
PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]

# Don't try to talk to a real SMTP server; collect mail in memory instead.
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"
