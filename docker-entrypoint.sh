#!/bin/bash

echo "+> pip list..."
pip list

CONTAINER_ROLE=${CONTAINER_ROLE:-"stack"}
echo "+> CONTAINER_ROLE is ${CONTAINER_ROLE}"

if [ "${CONTAINER_ROLE,,}" = "worker" ]; then

  # WORKER (at least one Pod will be a worker)
  echo "+> Running as WORKER"
  /code/launch-worker.sh

elif [ "${CONTAINER_ROLE,,}" = "beat" ]; then

  # BEAT (only one Pod will be a beat Pod)
  echo "+> Running as BEAT"
  /code/launch-beat.sh

elif [ "${CONTAINER_ROLE,,}" = "stack" ]; then

  # STACK (at least one Pod will be a stack)
  echo "+> Running as STACK"
  /code/launch-stack.sh

else
  
  # Unsupported Role! (this should never happen)
  echo "+> CONTAINER_ROLE is not supported!"

fi

echo "-> Exiting..."
