#!/bin/bash
echo "Running migrations..."
cd /code
# initial migrations for existing stack
python manage.py migrate auth
python manage.py migrate viewer
python manage.py migrate scoring hypothesis car

# make and apply new migrations not already in stack
python /code/manage.py makemigrations auth
python /code/manage.py migrate

python /code/manage.py makemigrations
python /code/manage.py migrate  # Apply database migrations - weird order is due to https://stackoverflow.com/questions/31417470/django-db-utils-programmingerror-relation-app-user-does-not-exist-during-ma

echo "Running collectstatic..."
python manage.py collectstatic --noinput -v 0 # collect static files
