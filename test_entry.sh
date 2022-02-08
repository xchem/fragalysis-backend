#!/bin/bash
/bin/bash /code/makemigrations.sh
export ISPYB_FLAG=""

python /code/manage.py test --no-input