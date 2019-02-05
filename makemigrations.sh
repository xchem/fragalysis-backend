#!/bin/bash
echo "Running migrations..."
cd /code
python /code/manage.py makemigrations auth
python /code/manage.py migrate
python /code/manage.py makemigrations xchem_db
python /code/manage.py migrate --fake-initial
python /code/manage.py makemigrations viewer
python /code/manage.py migrate  --fake-initial
python /code/manage.py makemigrations scoring hypothesis hotspots
python /code/manage.py migrate  --fake-initial
python /code/manage.py makemigrations
python /code/manage.py migrate  # Apply database migrations - weird order is due to https://stackoverflow.com/questions/31417470/django-db-utils-programmingerror-relation-app-user-does-not-exist-during-ma
echo "Running collectstatic..."
python manage.py collectstatic --noinput -v 0 # collect static files
