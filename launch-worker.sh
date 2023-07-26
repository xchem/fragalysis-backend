#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

echo "Running celery..."

CONCURRENCY=${WORKER_CONCURRENCY:-4}

export C_FORCE_ROOT=true
celery --app fragalysis worker \
    --concurrency ${CONCURRENCY}
