#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

# The tasks' scratch space (it is not shared: tasks pass storage names
# to each other, never local paths - see viewer/media_storage.py).
tmpdir="${TMPDIR:-/tmp}"
echo "Preparing tmp ($tmpdir)..."
mkdir -p ${tmpdir}

# part of debugging issue 1609, missing template protein.
# The watcher watches MEDIA_ROOT, so it only runs with the filesystem
# media storage (MEDIA_STORAGE, see fragalysis/settings.py and issue 1038).
media_storage="${MEDIA_STORAGE:-FILESYSTEM}"
if [ "${media_storage^^}" == "FILESYSTEM" ]; then
    echo "Starting media deletion watcher..."
    /code/filewatcher.sh &
else
    echo "Not starting the media deletion watcher (MEDIA_STORAGE=${MEDIA_STORAGE})"
fi

CONCURRENCY=${WORKER_CONCURRENCY:-4}

echo "Running celery (CONCURRENCY=${CONCURRENCY})..."
export C_FORCE_ROOT=true
celery --app fragalysis worker \
    --concurrency ${CONCURRENCY}
