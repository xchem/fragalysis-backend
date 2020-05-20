#!/bin/bash
/code/launch-stack.sh
echo "starting redis server..."
redis-server &
echo "starting celery daemon..."
celery -A fragalysis worker -l info &