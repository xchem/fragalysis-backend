#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

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

# nginx needs a DNS server to find the S3 host of the media storage (issue 1038).
# Use NGINX_RESOLVER (one or more addresses) if it is set,
# otherwise use the first nameserver in /etc/resolv.conf.
# The nginx resolver does not use search domains, so in Kubernetes
# the MEDIA_S3_ENDPOINT_URL host must be a fully qualified domain name.
if [ -n "${NGINX_RESOLVER:-}" ]; then
    resolver="${NGINX_RESOLVER}"
else
    resolver=$(awk '/^nameserver/ {print $2; exit}' /etc/resolv.conf)
    # nginx needs an IPv6 address in brackets
    if [[ "${resolver}" == *:* ]]; then
        resolver="[${resolver}]"
    fi
fi
if [ -n "${resolver}" ]; then
    echo "Using DNS server ${resolver} for nginx"
    echo "resolver ${resolver} valid=30s;" > /etc/nginx/frag_resolver.conf
else
    echo "WARNING: No DNS server for nginx (S3 media files will not be served)"
fi

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
