[![build dev](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-dev.yaml/badge.svg)](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-dev.yaml)
[![build staging](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-staging.yaml/badge.svg)](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-staging.yaml)
[![build production](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-production.yaml/badge.svg)](https://github.com/alanbchristie/fragalysis-backend/actions/workflows/build-production.yaml)

[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.0.1-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend/blob/master/LICENSE.txt)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis-backend.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis-backend/context:python)

# Fragalysis backend
Django server for Fragalysis with DRF API and loaders for data. Has components to serve React web-app.

## Documentation: https://fragalysis-backend.readthedocs.io/en/latest/index.html ##

# Dev environment setup
### Background

The stack consists of three services, running as containers: -

- a Postgres database (note: this used to be MySQL)
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

## Setup

**1.Create project directory, e.g.**


```
mkdir fragalysis
```

**2.Clone repositories inside your project's directory**


You can clone original `xchem` repositories or your forked e.g. `m2ms` and checkout any branch if necessary
```
git clone https://github.com/xchem/fragalysis-backend.git
git clone https://github.com/xchem/fragalysis-frontend.git
git clone https://github.com/InformaticsMatters/dls-fragalysis-stack-openshift.git
```
Note: 
To successful build, it should exist following directories from repository cloning.
Frontend is also important, because DJANGO server will serve this directory!
```$xslt
fragalysis/fragalysis-frontend/
fragalysis/fragalysis-backend/
fragalysis/dls-fragalysis-stack-openshift/
```

**3.Create some key data directories**


In `fragalysis/` directory run script to create data directory structure
```
mkdir -p data/input/django_data/EXAMPLE
mkdir -p data/neo4j/data
mkdir -p data/neo4j/logs
mkdir -p data/stack/media
mkdir -p data/stack/logs
mkdir -p data/media/compound_sets
mkdir -p data/postgre/data
```
**4.Populating database** 

Copy to `fragalysis/data/input/django_data/EXAMPLE` your PDB data, before you can launch the application.

If not exists file `fragalysis/data/input/django_data/EXAMPLE/TARGET_LIST` create it and content with list of your data, for example:
```
Mpro, NUDT7A,...
```

## Making database migrations
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

    pip install --upgrade pip
    pip install -r build-requirements.txt
    pre-commit install -t commit-msg -t pre-commit

Now the project's rules will run on every commit.

## Start
Start `Fragalysis stack` (All infrastructure - databases + populating data)

```
docker-compose up -d
```

`Please wait, it takes a minute until all containers are fully started.`

Test if we are running at [http://localhost:8080](http://localhost:8080)

If needed stop containers

```
docker-compose down
```


Note: The first run will be probably not successful. To fix after first run: 
- Delete the 'mysql_data' folder located in the data folder (Step 3) 
- Delete the migrations by running the commands below in the project directory (Step 1)
```
find . -path "*/migrations/*.py" -not -name "__init__.py" -delete
find . -path "*/migrations/*.pyc"  -delete
```
- In the 'fragalysis-frontend' folder rerun the docker compose command another two times 
```
docker-compose -f docker-compose.dev.yml up
```

### Populate database with CompoundSet
The first you need to have upload key.
#### Generate upload key
connect to web_dock service and run following script
```
python manage.py shell
from viewer.models import CSetKeys
new_key = CSetKeys()
new_key.name = 'test'
new_key.save()
print(new_key.uuid)
```

It is important to have key in following format (with `-`)
`f8c4ea0f-6b81-46d0-b85a-3601135756bc` 

#### Upload page
Visit `localhost:8080/viewer/upload_cset`

The target name is for example `Mpro`, select `sdf` file and you don't need a `pdb` file. 
Before upload generate upload key.

## Settings

### Development mode in DJANGO

In `settings.py`, the value of `DEBUG` is set to `False` by default. This can be set locally by adding the 
following line to the `web_dock` section of your docker_compose file:

```
    DEBUG_FRAGALYSIS: True
```

or adding the following line to the AWX template: 
```
    stack_debug_fragalysis: True
```
Note that this will also enable logging on the stack, so when coding you can use the logger and then jump on the 
webdock and check the logfile on /logs/ (although this doesn't work from celery yet).

### Sentry error logging

In `settings.py`, this is controlled by setting the value of `SENTRY_DNS`. To enable it, you need to set it to a valid
Sentry DNS entry. For Diamond, this can be set locally by adding the 
following line to the `web_dock` section of your docker_compose file:

```
    SENTRY_DNS: https://27fa0675f555431aa02ca552e93d8cfb@o194333.ingest.sentry.io/1298290
```

or adding the following line to the AWX template: 

```
    stack_backend_sentry_dsn: https://27fa0675f555431aa02ca552e93d8cfb@o194333.ingest.sentry.io/1298290
```

### Design Documents

As the application evolves several design documents have been written detailing improvements. These may be useful
for background reading on why decisions have been made.

The documents will be stored in the /design_docs folder in the repo. Current docs are listed below:
- [Fragalysis Discourse Design](design_docs/Fragalysis_Discourse_v0.2.pdf)
- [Fragalysis Tags Design V1.0](design_docs/Fragalysis_Tags_Design_V1.0.pdf)
- [Fragalysis Design #651 Fix Data Download V2.0](design_docs/Fragalysis_Design_651_Fix_Data_Download_V2.0.pdf)
- [Fragalysis Job Launcher V1.0](design_docs/Fragalysis_Job_Launcher_V1.0.pdf)
- [Fragalysis Job Launcher V2.0](design_docs/Fragalysis_Job_Launcher_Phase2_V1.0.pdf)

[pre-commit]: https://pre-commit.com
