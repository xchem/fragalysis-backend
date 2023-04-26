# Fragalysis backend

![GitHub release (latest by date)](https://img.shields.io/github/v/release/xchem/fragalysis-backend)

[![build dev](https://github.com/xchem/fragalysis-backend/actions/workflows/build-dev.yaml/badge.svg)](https://github.com/xchem/fragalysis-backend/actions/workflows/build-dev.yaml)
[![build staging](https://github.com/xchem/fragalysis-backend/actions/workflows/build-staging.yaml/badge.svg)](https://github.com/xchem/fragalysis-backend/actions/workflows/build-staging.yaml)
[![build production](https://github.com/xchem/fragalysis-backend/actions/workflows/build-production.yaml/badge.svg)](https://github.com/xchem/fragalysis-backend/actions/workflows/build-production.yaml)

[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend/blob/master/LICENSE.txt)

The Django server for Fragalysis with Django REST Framework for the API
and loaders for data.

>   See additional documentation relating to the backend on ReadTheDocs at
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
- Docker Compose
- Git

## Local development

>   These _local development_ notes are **Deprecated**, please see the
    documentation relating to the [Kubernetes Stack] deployment on ReadTheDocs where
    the development architecture is described in detail.

Create project directory: -

    mkdir fragalysis

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

Start the stack and generate an upload key: -

    docker-compose up -d

>   You may want to use a `.env` file to set some of the environment values defined in the
    `docker-compose.yml` file. The `.env` is ignored by git but allows you to set a number
    of variables, some of which may be _sensitive_.

Connect to `stack` service and run following script: -

    python manage.py shell
    from viewer.models import CSetKeys
    new_key = CSetKeys()
    new_key.name = 'test'
    new_key.save()
    print(new_key.uuid)

It is important to have key in following format: `f8c4ea0f-6b81-46d0-b85a-3601135756bc` 

With the stack running, upload a compound set by visiting `localhost:8080/viewer/upload_cset`

The target name is for example `Mpro`, select `sdf` file and you don't need a `pdb` file. 

## Database migrations
The best approach is to spin-up the development stack (locally) using
`docker-compose` and then shell into the Django (stack). For example,
to make new migrations called "add_job_request_start_and_finish_times"
for the viewer's models run the following: -

>   Before starting postgres, if you need to, remove any pre-existing db (if one exists)
    with `rm -rf ../data` (a directory maintained above the stack clone)

    docker-compose up -d
    docker-compose exec stack bash

Then from within the stack make the migrations. Here we're migrating the `viewer`
application...

    python manage.py makemigrations viewer --name "add_job_request_start_and_finish_times"

Exit the container and tear-down the deployment: -

    docker-compose down

>   The migrations will be written to your clone's filesystem as the clone directory
    is mapped into the container as a volume. You just need to commit the
    migrations to Git.

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

[kubernetes stack]: https://dls-fragalysis-stack-kubernetes.readthedocs.io/en/latest/index.html#
[pre-commit]: https://pre-commit.com
[readthedocs]: https://fragalysis-backend.readthedocs.io/en/latest/index.html
