#!/bin/bash

echo "+> CONTAINER_ROLE is ${CONTAINER_ROLE}"

if [ "${CONTAINER_ROLE,,}" = "worker" ]; then

  # WORKER
  #
  # This container is expected to act a
  # Fragalysis Stack Celery Task Worker
  echo "+> Running as WORKER"
  /code/launch-worker.sh

else

  # STACK
  #
  # This image is expected to act as the
  # Fragalysis Stack _first responder_ - i.e. the REST endpoint.
  echo "+> Running as STACK"
  /code/launch-stack.sh

fi
