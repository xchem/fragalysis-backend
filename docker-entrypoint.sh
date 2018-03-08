#!/bin/bash
npm install
python manage.py makemigrations auth
python manage.py migrate auth
python manage.py makemigrations scoring
python manage.py migrate scoring
python manage.py makemigrations
python manage.py migrate     # Apply database migrations - weird order is due to https://stackoverflow.com/questions/31417470/django-db-utils-programmingerror-relation-app-user-does-not-exist-during-ma
python loader.py
python manage.py collectstatic --clear --noinput -v 0 # clearstatic files
python manage.py collectstatic --noinput -v 0 # collect static files
# Prepare log files and start outputting logs to stdout
touch /srv/logs/gunicorn.log
touch /srv/logs/access.log
tail -n 0 -f /srv/logs/*.log &
# Start Gunicorn processes
echo Starting Gunicorn.
exec gunicorn fragalysis.wsgi:application \
    --name fragalysis \
    --bind unix:django_app.sock \
    --workers 3 \
    --log-level=info \
    --log-file=/srv/logs/gunicorn.log \
    --access-logfile=/srv/logs/access.log &
exec service nginx start
