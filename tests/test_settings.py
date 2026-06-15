# pylint: disable=wildcard-import,unused-wildcard-import
"""Django settings used when running the test suite under pytest.

This module imports everything from the normal application settings and then
applies test-only overrides. It is selected via ``DJANGO_SETTINGS_MODULE`` in
``[tool.pytest.ini_options]`` (see ``pyproject.toml``).

The schema requires a real PostgreSQL instance (pgvector, simple_history, etc.)
- SQLite cannot stand in for it - but it does not require the backend container
image. By default we connect to a PostgreSQL on ``localhost`` (the ``database``
service of ``docker-compose.test.yml``, or a CI service container), so the suite
runs with a bare ``pytest`` on the host. The ``POSTGRESQL_*`` environment
variables still override these defaults when supplied.
"""
import os
import tempfile

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
# (Project.user_id membership). TA_AUTH_SERVICE is inherited from the real
# settings (empty unless the env var is set) so these tests also guard against
# that setting going missing. Tests that want the external TA-authenticator
# path set settings.TA_AUTH_SERVICE themselves (e.g. via the pytest-django
# `settings` fixture) and patch ta_auth_connector.get_auth_target_access.

# Use the plain PostgreSQL backend for tests instead of the django_prometheus
# wrapper. The wrapper registers process-wide Prometheus collectors at import
# time, which clashes ("Duplicated timeseries in CollectorRegistry") when
# settings are imported repeatedly under pytest. The schema is identical.
DATABASES["default"]["ENGINE"] = "django.db.backends.postgresql"  # noqa: F405

# Connect to a localhost PostgreSQL by default (the application settings default
# to the docker-network name "database"/user "fragalysis"). This lets the suite
# run on the host against the docker-compose.test.yml "database" container - or
# a CI service container - without the backend image. An explicit POSTGRESQL_*
# environment variable still wins.
DATABASES["default"]["HOST"] = os.environ.get(
    "POSTGRESQL_HOST", "127.0.0.1"
)  # noqa: F405
DATABASES["default"]["USER"] = os.environ.get(
    "POSTGRESQL_USER", "postgres"
)  # noqa: F405

# Fast, insecure password hashing - we never need real hashing strength in
# tests and this noticeably speeds up user creation.
PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]

# Don't try to talk to a real SMTP server; collect mail in memory instead.
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"

# Tests that touch the media tree (create/delete sub-directories) need a
# writable MEDIA_ROOT. The application default ("/code/media/") only exists in
# the container, so use an isolated temp directory that also works on the host.
MEDIA_ROOT = tempfile.mkdtemp(prefix="fragalysis-test-media-")

# The application logging config writes to rotating files under BASE_DIR/logs,
# a directory that only exists in the container (a mounted volume). Redirect
# those files to a writable temp directory so the suite runs on a bare host / CI
# checkout. (Skipped when DISABLE_LOGGING_FRAMEWORK leaves LOGGING undefined.)
if "LOGGING" in globals():
    _test_log_dir = tempfile.mkdtemp(prefix="fragalysis-test-logs-")
    _handlers: dict = LOGGING["handlers"]  # type: ignore[assignment]  # noqa: F405
    for _handler in _handlers.values():
        if "filename" in _handler:
            _handler["filename"] = os.path.join(
                _test_log_dir, os.path.basename(_handler["filename"])
            )
