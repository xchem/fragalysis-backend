#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

# part of debugging issue 1609, missing template protein
if [ "${WATCH_MEDIA_DELETIONS:-false}" = "true" ]; then
    echo "Starting media deletion watcher..."
    /code/filewatcher.sh &
fi

CONCURRENCY=${WORKER_CONCURRENCY:-4}

echo "Running celery (CONCURRENCY=${CONCURRENCY})..."
export C_FORCE_ROOT=true
celery --app fragalysis worker \
    --concurrency ${CONCURRENCY}
