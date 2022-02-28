#!/bin/bash
/bin/bash /code/makemigrations.sh
export ISPYB_FLAG=""

python3 /code/manage.py test --no-input