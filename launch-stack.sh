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

echo "Starting Gunicorn..."
gunicorn --config gunicorn.conf.py fragalysis.wsgi:application

# added as a fix to #1215, mixing http and https requests to enable
# local development with http. Need to set the env variable in compose
# file or .env.

# NB! this is probably a workaround for some other issue we haven't
# discovered yet. It suddenly broke but right now seems to work in
# firefox. 12Hopefully won't be necessary soon
echo proxy_set_header X-Forwarded-Proto "${PROXY_FORWARDED_PROTO_HEADER:-https};"  >> /etc/nginx/frag_proxy_params

# Render nginx config templates. Default of 3600s (1 hour) supports
# large file uploads (#942); override with NGINX_TIMEOUT_S.
# Only ${NGINX_TIMEOUT_S} is substituted so nginx's own $vars
# (e.g. $http_host, regex anchors) pass through untouched.
export NGINX_TIMEOUT_S="${NGINX_TIMEOUT_S:-3600}"
echo "Rendering nginx config (NGINX_TIMEOUT_S=${NGINX_TIMEOUT_S})..."
envsubst '${NGINX_TIMEOUT_S}' \
    < /etc/nginx/templates/nginx.conf.template \
    > /etc/nginx/nginx.conf
envsubst '${NGINX_TIMEOUT_S}' \
    < /etc/nginx/templates/default.conf.template \
    > /etc/nginx/sites-available/default.conf

echo "Testing nginx config..."
nginx -tq

if [ "${HOSTNAME}" = "stack-0" ]; then
    echo "Clearing the page cache"
    python manage.py clear_cache
    echo "Launching service health check queries"
    python manage.py start_service_queries &
    echo "Launching download cleanup scheduler"
    python manage.py start_download_cleanup &
fi

echo "Running nginx..."
nginx
