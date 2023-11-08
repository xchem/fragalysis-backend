#!/usr/bin/env bash

# A convenient script to build and push the docker image to the registry.
# Used for lightweight local dev. CI uses its own process.
#
# Before running, set the following environment variables:
#
#   export BE_IMAGE_TAG=1187.1
#   export BE_NAMESPACE=alanbchristie
#
# Alan Christie
# Nov 2023

# Exit conditions...
# -x to diplay commands as they are run (expanding variables)
# -e exits on error
# -o (for option) pipefail exits on command pipe failures
set -xeo pipefail

: "${BE_IMAGE_TAG?Need to set BE_IMAGE_TAG}"
: "${BE_NAMESPACE?Need to set BE_NAMESPACE}"

poetry export --without-hashes --without dev --output requirements.txt
docker build --no-cache . -t ${BE_NAMESPACE}/fragalysis-backend:${BE_IMAGE_TAG}
docker push ${BE_NAMESPACE}/fragalysis-backend:${BE_IMAGE_TAG}
