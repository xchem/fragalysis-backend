#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

# part of debugging issue 1609, missing template protein
echo "Starting media deletion watcher..."
/code/filewatcher.sh &

echo "Running migrations..."
cd /code
python manage.py migrate

echo "Loading fixtures..."
python manage.py loaddata tagcategories.json

echo "Running collectstatic..."
python manage.py collectstatic --noinput -v 0

echo "Creating superuser..."
# Automatically create the superuser...
script="
from django.contrib.auth.models import User
username = '$WEB_DJANGO_SUPERUSER_NAME'
password = '$WEB_DJANGO_SUPERUSER_PASSWORD'
email = '$WEB_DJANGO_SUPERUSER_EMAIL'
if not username or not password:
    username = 'admin'
    password = 'UNSECURED'
if User.objects.filter(username=username).count()==0:
    User.objects.create_superuser(username, email, password)
    print(f'Superuser created ({username}|{password}).')
else:
    print('Superuser creation skipped - an admin user already exists.')
"
printf "$script" | python manage.py shell

echo "Preparing logging..."
touch /code/logs/gunicorn.log
touch /code/logs/access.log
touch /code/logs/backend.log

CONCURRENCY=${STACK_CONCURRENCY:-4}

echo "Starting Gunicorn (CONCURRENCY=${CONCURRENCY})..."
gunicorn fragalysis.wsgi:application \
    --daemon \
    --name fragalysis \
    --bind unix:django_app.sock \
    --timeout 3000 \
    --workers ${CONCURRENCY} \
    --log-level=debug \
    --log-file /code/logs/gunicorn.log \
    --error-logfile /code/logs/gunicorn-err.log \
    --access-logfile /code/logs/gunicorn-access.log

# added as a fix to #1215, mixing http and https requests to enable
# local development with http. Need to set the env variable in compose
# file or .env.

# NB! this is probably a workaround for some other issue we haven't
# discovered yet. It suddenly broke but right now seems to work in
# firefox. 12Hopefully won't be necessary soon
echo proxy_set_header X-Forwarded-Proto "${PROXY_FORWARDED_PROTO_HEADER:-https};"  >> /etc/nginx/frag_proxy_params

echo "Testing nginx config..."
nginx -tq

if [ "${HOSTNAME}" = "stack-0" ]; then
  echo "Launching service health check queries"
  python manage.py start_service_queries &
  echo "Launching download cleanup scheduler"
  python manage.py start_download_cleanup &
fi

echo "Running nginx..."
nginx
