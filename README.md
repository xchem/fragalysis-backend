[![build dev](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-dev.yaml/badge.svg)](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-dev.yaml)
[![build staging](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-staging.yaml/badge.svg)](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-staging.yaml)
[![build production](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-production.yaml/badge.svg)](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-production.yaml)

[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.0.1-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend/blob/master/LICENSE.txt)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis-backend.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis-backend/context:python)

# Fragalysis backend
Django server for Fragalysis with DRF API and loaders for data. Has components to serve React web-app.

>   See additional documentation relating to the Stack, especially the development
    process and deployment mechanism on ReadTheDocs at
    https://fragalysis-backend.readthedocs.io/en/latest/index.html

## Background

The stack consists of three services, running as containers: -

- a Postgres database
- a neo4j graph database
- the fraglaysis stack

The stack is formed from code resident in a number of repositories.
Begin by forking repositories you anticipate editing (although you really want
to consider forking all the repositories as this is a relatively low-cost
operation).

The repositories are:

- [xchem/fragalysis](https://github.com/xchem/fragalysis)
- [xchem/fragalysis-frontend](https://github.com/xchem/fragalysis-frontend)
- [xchem/fragalysis-backend](https://github.com/xchem/fragalysis-backend)
- [xchem/fragalysis-stack](https://github.com/xchem/fragalysis-stack)

### Prerequisites

- Docker
- Docker-compose
- Git

## Local development

Create project directory: -

mkdir fragalysis
mkdir fragalysis
```
    mkdir fragalysis
```


Clone repositories inside your project's directory: -


>   You can clone original `xchem` repositories or your forked
    e.g. `m2ms` and checkout any branch if necessary

    git clone https://github.com/xchem/fragalysis-backend.git
    git clone https://github.com/xchem/fragalysis-frontend.git
    git clone https://github.com/InformaticsMatters/dls-fragalysis-stack-openshift.git

Note: To successfully build and run, the following directories should exit
from repository cloning. The frontend is also important, because the django
server will serve this directory!

    fragalysis/fragalysis-frontend/
    fragalysis/fragalysis-backend/
    fragalysis/dls-fragalysis-stack-openshift/

Create some key data directories: -

In the `fragalysis/` directory run the following to create a data directory structure: -

    mkdir -p data/input/django_data/EXAMPLE
    mkdir -p data/neo4j/data
    mkdir -p data/neo4j/logs
    mkdir -p data/stack/media
    mkdir -p data/stack/logs
    mkdir -p data/media/compound_sets
    mkdir -p data/postgre/data

Populating database: -

Copy to `fragalysis/data/input/django_data/EXAMPLE` your PDB data, before you can launch the application.

If the file `fragalysis/data/input/django_data/EXAMPLE/TARGET_LIST` does not exist,
create it and add a content list of your data, for example:

    Mpro, NUDT7A,...

## Database migrations
The best approach is to spin-up the development stack (locally) using
`docker-compose` and then shell into the Django (stack). For example,
to make new migrations called "add_job_request_start_and_finish_times"
for the viewer's models run the following: -

>   Before starting postgres, if you need to, remove any pre-existing Db (if one exists)
    with `rm -rf data`

    docker-compose up -d
    docker-compose exec stack bash

Then from within the stack make the migrations. Here we're migrating the `viewer`
application...

    python manage.py makemigrations viewer --name "add_job_request_start_and_finish_times"

Exit the container and tear-down the deployment: -

    docker-compose down

## Pre-commit hooks
The project uses [pre-commit] to enforce linting of files prior to committing
them to the upstream repository.

To get started review the pre-commit utility and then set-up your local clone
by following the **Installation** and **Quick Start** sections of the
pre-commit documentation.

Ideally from a Python environment...

    python -m venv venv
    source venv/bin/activate

    pip install --upgrade pip
    pip install -r build-requirements.txt
    pre-commit install -t commit-msg -t pre-commit

Now the project's rules will run on every commit and you can check the
state of the repository as it stands with...

    pre-commit run --all-files

## Start
Start `Fragalysis stack` (All infrastructure - databases + populating data)

### Generate upload key
Connect to `stack` service and run following script

    python manage.py shell
    from viewer.models import CSetKeys
    new_key = CSetKeys()
    new_key.name = 'test'
    new_key.save()
    print(new_key.uuid)

It is important to have key in following format (with `-`)
`f8c4ea0f-6b81-46d0-b85a-3601135756bc` 

#### Upload page
Visit `localhost:8080/viewer/upload_cset`

The target name is for example `Mpro`, select `sdf` file and you don't need a `pdb` file. 
Before upload generate upload key.

## Sentry error logging

In `settings.py`, this is controlled by setting the value of `SENTRY_DNS`.
To enable it, you need to set it to a valid Sentry DNS entry.
For Diamond, this can be set locally by adding the 
following line to the `stack` section of your docker_compose file:

    SENTRY_DNS: https://<SENTRY_DNS>

## Design Documents

As the application evolves several design documents have been written detailing
improvements. These may be useful for background reading on why decisions have been made.

The documents will be stored in the /design_docs folder in the repo. Current docs are listed below: -

- [Fragalysis Discourse Design](design_docs/Fragalysis_Discourse_v0.2.pdf)
- [Fragalysis Tags Design V1.0](design_docs/Fragalysis_Tags_Design_V1.0.pdf)
- [Fragalysis Design #651 Fix Data Download V2.0](design_docs/Fragalysis_Design_651_Fix_Data_Download_V2.0.pdf)
- [Fragalysis Job Launcher V1.0](design_docs/Fragalysis_Job_Launcher_V1.0.pdf)
- [Fragalysis Job Launcher V2.0](design_docs/Fragalysis_Job_Launcher_Phase2_V1.0.pdf)

---

[pre-commit]: https://pre-commit.com
