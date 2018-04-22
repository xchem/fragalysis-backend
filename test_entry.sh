#!/bin/bash
python /code/manage.py makemigrations auth
python /code/manage.py migrate auth
python /code/manage.py makemigrations scoring pandda hypothesis
python /code/manage.py migrate scoring pandda hypothesis
python /code/manage.py makemigrations
python /code/manage.py migrate
python /code/manage.py test --no-input