#!/bin/bash
export ISPYB_FLAG=""

# Run tests that don't require a 'test' database - i.e. unit-level tests
python manage.py test --settings tests.nodb_settings --tag nodb
# Test everything else...
python manage.py test --exclude-tag nodb --no-input
