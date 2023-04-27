# The Fragalysis Stack (Backend)

![GitHub release (latest by date)](https://img.shields.io/github/v/release/xchem/fragalysis-backend)

[![build dev](https://github.com/xchem/fragalysis-backend/actions/workflows/build-dev.yaml/badge.svg)](https://github.com/xchem/fragalysis-backend/actions/workflows/build-dev.yaml)
[![build staging](https://github.com/xchem/fragalysis-backend/actions/workflows/build-staging.yaml/badge.svg)](https://github.com/xchem/fragalysis-backend/actions/workflows/build-staging.yaml)
[![build production](https://github.com/xchem/fragalysis-backend/actions/workflows/build-production.yaml/badge.svg)](https://github.com/xchem/fragalysis-backend/actions/workflows/build-production.yaml)

[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend/blob/master/LICENSE.txt)

The Django server for Fragalysis using the Django REST Framework ([DRF])
for the API, and loaders for data.

>   See additional documentation relating to the backend on ReadTheDocs at
    https://fragalysis-backend.readthedocs.io/en/latest/index.html

## Background
The stack consists of three services, running as containers: -

- a Postgres database
- a neo4j graph database
- the Fraglaysis "stack"

The stack is formed from code resident in a number of repositories.
Begin by forking repositories you anticipate editing (although you really want
to consider forking all the repositories as this is a relatively low-cost
operation).

The repositories are: -

- [xchem/fragalysis](https://github.com/xchem/fragalysis)
- [xchem/fragalysis-frontend](https://github.com/xchem/fragalysis-frontend)
- [xchem/fragalysis-backend](https://github.com/xchem/fragalysis-backend)
- [xchem/fragalysis-stack](https://github.com/xchem/fragalysis-stack)

The stack is deployed as a container image to [Kubernetes] using [Ansible] playbooks
that can be found in the Ansible repository, where additional development and deployment
documentation can also be found: -

- [informaticsmatters/dls-fragalysis-stack-kubernetes](https://github.com/InformaticsMatters/dls-fragalysis-stack-kubernetes)

## Command-line access to the API
With the stack (or backend) running you should be able to access the REST API. From
the command-line you can use curl or [httpie]. Here's we use the httpie utility to
`GET` the API root (which does not require authentication)...

    http :8080/api/

To use much of the remainder of the API you will need to authenticate. but some endpoints allow you to use a token,
obtained from the corresponding Keycloak authentication service. If you are
running a local stack a client ID exists that should work for you, assuming you have
a Keycloak user identity. So, with a few variables: -

    TOKEN_URL=keycloak.example.com/auth/realms/xchem/protocol/openid-connect/token
    CLIENT_ID=fragalysis-local
    CLIENT_SECRET=00000000-0000-0000-0000-000000000000
    USER=someone
    PASSWORD=password123

...you can then get an API token. Here we're using `httpie`and `jq`: -

    TOKEN=$(http --form POST https://$TOKEN_URL/ \
        grant_type=password \
        client_id=$CLIENT_ID \
        client_secret=$CLIENT_SECRET \
        username=$USER \
        password=$PASSWORD | jq -r '.access_token')

The token should last for at least 15 minutes, depending on the Keycloak configuration.
With the Token you should then be able to make authenticated requests to the API on your
local stack: -

    ENDPOINT=api/compound-identifier-types
    http :8080/$ENDPOINT/ "Authorization:Bearer $TOKEN"
    
## Logging
The backend writes log information in the container to `/code/logs/backend.log`. This is
typically persisted between container restarts on Kubernetes with a separate volume mounted
at `/code/logs`.

## Database migrations
The best approach is to spin-up the development stack (locally) using
`docker-compose` and then shell into the Django (stack). For example,
to make new migrations called "add_job_request_start_and_finish_times"
for the viewer's models run the following: -

>   Before starting postgres, if you need to, remove any pre-existing database (if one exists)
    with `rm -rf ../data` (a directory maintained above the repository's clone)

    docker-compose up -d

Then from within the stack container make the migrations. Here we're migrating the `viewer`
application...

    docker-compose exec stack bash

    python manage.py makemigrations viewer --name "add_job_request_start_and_finish_times"

Exit the container and tear-down the deployment: -

    docker-compose down

>   The migrations will be written to your clone's filesystem as the clone directory
    is mapped into the container as a volume. You just need to commit the
    migrations that have been written to the local directory to Git.

## Sentry error logging
In `settings.py`, this is controlled by setting the value of `SENTRY_DNS`.
To enable it, you need to set it to a valid Sentry DNS entry.
For Diamond, this can be set locally by adding the 
following line to the `stack` section of your docker_compose file:

    SENTRY_DNS: https://<SENTRY_DNS>

## Compiling the documentation
Because the documentation uses Sphinx and its `autodoc` module, compiling the
documentation needs all the application requirements. As this is often impractical
on the command-line, the most efficient way to build the documentation is from within
the stack container (as shown below).

    docker-compose up -d
    docker-compose exec stack bash

    pip install sphinx==5.3.0
    pip install importlib-metadata~=4.0
    
    cd docs
    sphinx-build -b html source/ build/

>   The current version of Python used in the Django container image is **3.7**
    and this suffers from an import error relating to celery. It is fixed by
    using a pre-v5.0 version of `importlib-metadata` as illustrated in the above example.
    (see https://stackoverflow.com/questions/73933432/)

The code directory is mounted in the container the documentation can then be committed
from the host machine.

## Local development
>   These _local development_ notes are **deprecated**, please see the
    documentation relating to the [Kubernetes Stack] deployment on ReadTheDocs where
    the development architecture is described in detail.

Prerequisites: -

- Docker
- Docker Compose
- Git
- Ideally a Linux host (although Windows and Mac should work)

Other 'handy' tools to have while developing are: -

- [jq]
- [httpie]

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

Build and start the stack and generate an upload key: -

>   You may want to use a `.env` file to set some of the environment values defined in the
    `docker-compose.yml` file. The `.env` is ignored by git but allows you to set a number
    of variables, some of which may be _sensitive_.

To run successfully you wil need to provide some variables, ideally in a `.env` file.
Namely: -

    OIDC_KEYCLOAK_REALM=https://keycloak.example.com/auth/realms/xchem
    OIDC_RP_CLIENT_ID=fragalysis-local
    OIDC_RP_CLIENT_SECRET=<secret>

Then run: -

    docker-compose build
    docker-compose up -d

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

## Design Documents
As the application has evolved several design documents have been written detailing
improvements. These may be useful for background reading on why decisions have been made.

The documents will be stored in the `/design_docs` folder in the repo.
Current docs are listed below: -

- [Fragalysis Discourse Design](design_docs/Fragalysis_Discourse_v0.2.pdf)
- [Fragalysis Tags Design V1.0](design_docs/Fragalysis_Tags_Design_V1.0.pdf)
- [Fragalysis Design #651 Fix Data Download V2.0](design_docs/Fragalysis_Design_651_Fix_Data_Download_V2.0.pdf)
- [Fragalysis Job Launcher V1.0](design_docs/Fragalysis_Job_Launcher_V1.0.pdf)
- [Fragalysis Job Launcher V2.0](design_docs/Fragalysis_Job_Launcher_Phase2_V1.0.pdf)

## Pre-commit
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

---

[ansible]: https://github.com/ansible/ansible
[drf]: https://www.django-rest-framework.org
[httpie]: https://pypi.org/project/httpie/
[jq]: https://stedolan.github.io/jq/
[kubernetes]: https://kubernetes.io
[kubernetes stack]: https://dls-fragalysis-stack-kubernetes.readthedocs.io/en/latest/index.html#
[pre-commit]: https://pre-commit.com
[readthedocs]: https://fragalysis-backend.readthedocs.io/en/latest/index.html
