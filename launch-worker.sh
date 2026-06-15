#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

tmpdir="${TMPDIR:-/code/media/tmp}"
echo "Preparing tmp ($tmpdir)..."
mkdir -p ${tmpdir}

# part of debugging issue 1609, missing template protein
echo "Starting media deletion watcher..."
/code/filewatcher.sh &

CONCURRENCY=${WORKER_CONCURRENCY:-4}

echo "Running celery (CONCURRENCY=${CONCURRENCY})..."
export C_FORCE_ROOT=true
celery --app fragalysis worker \
    --concurrency ${CONCURRENCY}
