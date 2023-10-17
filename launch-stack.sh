#!/bin/bash

# Exit conditions...
# -e exits on error,
# -o (for option) pipefail exits on command pipe failures
set -eo pipefail

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
    print('Superuser created.')
else:
    print('Superuser creation skipped.')
"
printf "$script" | python manage.py shell

echo "Preparing logging..."
touch /srv/logs/gunicorn.log
touch /srv/logs/access.log
touch /code/logs/logfile.log

echo "Starting Gunicorn..."
gunicorn fragalysis.wsgi:application \
    --daemon \
    --name fragalysis \
    --bind unix:django_app.sock \
    --timeout 3000 \
    --workers 3 \
    --log-level=debug \
    --log-file=/srv/logs/gunicorn.log \
    --access-logfile=/srv/logs/access.log

echo "Testing nginx config..."
nginx -t

echo "Running nginx..."
nginx
