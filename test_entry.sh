#!/bin/bash
/bin/bash /code/makemigrations.sh
export ISPYB_FLAG=""
# Dummy env variables for test
export IBM_API_KEY=123
export MANIFOLD_API_KEY=123
export MCULE_API_KEY=123
export SENDGRID_API_KEY=123

python /code/manage.py test --no-input