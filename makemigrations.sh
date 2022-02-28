#!/bin/bash
echo "Running migrations..."
cd /code
# initial migrations for existing stack
python3 manage.py migrate auth
python3 manage.py migrate viewer
python3 manage.py migrate scoring hypothesis car

# make and apply new migrations not already in stack
python3 /code/manage.py makemigrations auth
python3 /code/manage.py migrate

python3 /code/manage.py makemigrations
python3 /code/manage.py migrate  # Apply database migrations - weird order is due to https://stackoverflow.com/questions/31417470/django-db-utils-programmingerror-relation-app-user-does-not-exist-during-ma

echo "Running collectstatic..."
python3 manage.py collectstatic --noinput -v 0 # collect static files
