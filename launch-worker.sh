#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

CONCURRENCY=${WORKER_CONCURRENCY:-4}

echo "Running celery (CONCURRENCY=${CONCURRENCY})..."
export C_FORCE_ROOT=true
celery --app fragalysis worker \
    --concurrency ${CONCURRENCY}
