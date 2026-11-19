#!/bin/bash

# The image has no pip, so list the installed packages with importlib.metadata
echo "+> Installed packages..."
python -c 'from importlib.metadata import distributions
for name, version in sorted((d.name, d.version) for d in distributions()):
    print(name, version)'

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
