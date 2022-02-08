#!/bin/bash
/bin/bash /code/makemigrations.sh
export ISPYB_FLAG=""
# Dummy env variables for test
export IBM_API_KEY=""
export MANIFOLD_API_KEY=""
export MCULE_API_KEY=""
export SENDGRID_API_KEY=""

python /code/manage.py test --no-input