#!/bin/bash
echo "Running migrations..."
cd /code
# initial migrations for existing stack
python manage.py migrate viewer 0001_initial
python manage.py migrate viewer 0002_auto_20200615_1322
python manage.py migrate auth
python manage.py migrate xchem_db
python manage.py migrate scoring hypothesis hotspots
# make and apply new migrations not already in stack
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
