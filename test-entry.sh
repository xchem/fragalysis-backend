#!/bin/bash
# Entry point used by docker-compose.test.yml (and CI) to run the test suite.
# pytest configuration lives in pyproject.toml ([tool.pytest.ini_options]),
# including the test settings module and coverage options.
set -e

cd /code
pytest
