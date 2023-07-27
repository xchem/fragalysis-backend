#!/bin/bash

CONTAINER_ROLE=${CONTAINER_ROLE:-"stack"}
echo "+> CONTAINER_ROLE is ${CONTAINER_ROLE}"

if [ "${CONTAINER_ROLE,,}" = "worker" ]; then

  # WORKER (at least one Pod will be a worker)
  #
  # This container is expected to act a
  # Fragalysis Stack Celery Task Worker
  echo "+> Running as WORKER"
  /code/launch-worker.sh

elif [ "${CONTAINER_ROLE,,}" = "beat" ]; then

  # BEAT (only one Pod will be a beat Pod)
  #
  # This container is expected to act a
  # Fragalysis Stack Celery Beat Scheduler
  echo "+> Running as BEAT"
  /code/launch-beat.sh

else

  # STACK (at least one Pod will be a stack)
  #
  # This image is expected to act as the
  # Fragalysis Stack _first responder_ - i.e. the REST endpoint.
  echo "+> Running as STACK"
  /code/launch-stack.sh

fi
