#!/bin/bash
# Run the backend unit tests on the host with pytest - no backend image build.
#
# The test schema needs a real PostgreSQL (pgvector, simple_history), so this
# starts just the database container (not the backend image) and points pytest
# at it. Any arguments are passed straight through to pytest, for example:
#
#   ./run-unit-tests.sh
#   ./run-unit-tests.sh -k expired
#   ./run-unit-tests.sh viewer/tests/test_download_capacity.py
#
# First-time setup (once): install the host virtualenv with
#   poetry install --only main,test
set -euo pipefail

cd "$(dirname "$0")"

# Start only the PostgreSQL service (defined in docker-compose.test.yml) and
# wait until it reports healthy. This neither builds nor needs the backend image.
docker compose -f docker-compose.test.yml up --detach --wait database

# Ensure pytest (the 'test' dependency group) is present in the host virtualenv.
if ! poetry run python -c "import pytest" >/dev/null 2>&1; then
    echo "Installing test dependencies into the Poetry virtualenv..."
    poetry install --only main,test --no-root --no-directory
fi

# Connection details for the localhost database container. These match the
# defaults in tests/test_settings.py and are set here only to be explicit.
export POSTGRESQL_HOST=127.0.0.1
export POSTGRESQL_PORT=5432
export POSTGRESQL_DATABASE=frag
export POSTGRESQL_USER=postgres
export POSTGRESQL_PASSWORD=fragalysis
export CACHE_ENABLED=No

poetry run pytest "$@"
