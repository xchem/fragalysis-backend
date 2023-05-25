# The Fragalysis Backend

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
The **Backend** is part of the **Stack**, which consists of three services: -

- a Postgres database
- a neo4j graph database
- the Fraglaysis "stack"

The stack is formed from code resident in a number of repositories.
This one, and: -

- [xchem/fragalysis-frontend](https://github.com/xchem/fragalysis-frontend)
- [xchem/fragalysis-stack](https://github.com/xchem/fragalysis-stack)

Other, significant, repositories include: -

- [xchem/fragalysis](https://github.com/xchem/fragalysis)
- [xchem/fragalysis-api](https://github.com/xchem/fragalysis-api)
  
The stack is deployed as a container images to [Kubernetes] using [Ansible] playbooks
that can be found in the Ansible repository. Additional development and deployment
documentation can be found in the [informaticsmatters/dls-fragalysis-stack-kubernetes](https://github.com/InformaticsMatters/dls-fragalysis-stack-kubernetes) repository.

## Building and running (local)
The backend is a Docker container image and can be build and deployed locally using
`docker-compose`.

    docker-compose build

To run the application (which wil include deployment of the postgres and neo4j databases)
run: -

    docker-compose up -d

The postgres database is persisted in the `data` directory, outside of the repo.

You may need to provide a number of environment variables
that are employed in the container image. Fragalysis configuration depends on
a large number of variables, and the defaults may not be suitable for your needs.

The typical pattern with `docker-compose`, is to provide these variables in the
`docker-compose.yml` file and adjust their values (especially the sensitive ones)
using a local `.env` file (see [environment variables]).

The backend API, for example, should be available on port 8080 of your host
at `http://localhost:8080/api/`.

You can visit the `/accounts/login` endpoint to login (assuming you have setup the
appropriate environment variables for the container). This generates errors relating
to the fact that the FE/Webpack can’t be found. This looks alarming but you are logged in.

>   The backend no longer writes `.pyc` files (the `Dockerfile` sets the environment
    variable `PYTHONDONTWRITEBYTECODE`). This, and the fact the backend code is mapped
    into the container, allows you to make "live" changes to the code
    on your host and see them reflected in the container app without having to rebuild
    or restart the stack container.

When you want to spin-down the deployment run: -

    docker-compose down

## Command-line access to the API
With the stack (or backend) running you should be able to access the REST API. From
the command-line you can use [curl] or [httpie]. Here, we use `http` to
**GET** a response from the API root (which does not require authentication)...

    http :8080/api/

The response should contain a list of endpoint names and URLs, something like this...

```
{
    "action-type": "http://localhost:8080/api/action-type/",
    "cmpdchoice": "http://localhost:8080/api/cmpdchoice/",
    "cmpdimg": "http://localhost:8080/api/cmpdimg/",
    [...]
    "vector3ds": "http://localhost:8080/api/vector3ds/",
    "vectors": "http://localhost:8080/api/vectors/",
    "viewscene": "http://localhost:8080/api/viewscene/"
}
```

To use much of the remainder of the API you will need to authenticate.
Some endpoints allow you to use a token, obtained from the corresponding Keycloak
authentication service. If you are running a local stack a client ID exists that
should work for you, assuming you have a Keycloak user identity.
With a few variables: -

    TOKEN_URL=keycloak.example.com/auth/realms/xchem/protocol/openid-connect/token
    CLIENT_ID=fragalysis-local
    CLIENT_SECRET=00000000-0000-0000-0000-000000000000
    USER=someone
    PASSWORD=password123

...you should eb able to obtain an API token. Here we're using `http`and `jq`: -

    TOKEN=$(http --form POST https://$TOKEN_URL/ \
        grant_type=password \
        client_id=$CLIENT_ID \
        client_secret=$CLIENT_SECRET \
        username=$USER \
        password=$PASSWORD | jq -r '.access_token')

The token should last for at least 15 minutes, depending on the Keycloak configuration.
With the Token you should then be able to make authenticated requests to the API on your
local stack.

Here's an illustration of how to use the API from the command-line by getting, adding,
and deleting a `CompoundIdentifierType`: -

    ENDPOINT=api/compound-identifier-types

    http :8080/$ENDPOINT/ "Authorization:Bearer $TOKEN"
    RID=$(http post :8080/$ENDPOINT/ "Authorization:Bearer $TOKEN" name="XT345632" | jq -r '.id')
    http delete :8080/$ENDPOINT/$RID/ "Authorization:Bearer $TOKEN"

## Logging
The backend writes log information in the container to `/code/logs/backend.log`. This is
typically persisted between container restarts on Kubernetes with a separate volume mounted
at `/code/logs`.

## Database migrations
The best approach is to spin-up the development stack (locally) using
`docker-compose` and then shell into Django. For example,
to make new migrations called "add_job_request_start_and_finish_times"
for the viewer's model run the following: -

>   Before starting postgres, if you need to, remove any pre-existing database (if one exists)
    with `rm -rf ../data` (a directory maintained above the repository's clone)

    docker-compose up -d

Then from within the stack container make the migrations
(in this case for the `viewer`)...

    docker-compose exec stack bash

    python manage.py makemigrations viewer --name "add_job_request_start_and_finish_times"

Exit the container and tear-down the deployment: -

    docker-compose down

>   The migrations will be written to your clone's filesystem as the clone directory
    is mapped into the container as a volume. You just need to commit the
    migrations that have been written to the local directory to Git.

## Sentry error logging
[Sentry] can be used to log errors in the stack container image.

In `settings.py`, this is controlled by setting the value of `FRAGALYSIS_BACKEND_SENTRY_DNS`,
which is also exposed in the developer docker-compose file.
To enable it, you need to set it to a valid Sentry DNS value.

## Compiling the documentation
Because the documentation uses Sphinx and its `autodoc` module, compiling the
documentation needs all the application requirements. As this is often impractical
on the command-line, the most efficient way to build the documentation is from within
the stack container: -

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

The code directory is mounted in the container so the compiled documentation
can then be committed from the host machine.

## Pre-commit
The project uses [pre-commit] to enforce linting of files prior to committing
them to the upstream repository.

>   As fragalysis is a complex code-based (that's been maintained by a number of
    key developers) we currently limit the linting to the `viewer` application
    (see the `.pre-commit-config.yaml` file for details).
    In the future we might extend this to the entire code-base.

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

## Design documents
As the application has evolved several design documents have been written detailing
improvements. These may be useful for background reading on why decisions have been made.

The documents will be stored in the `/design_docs` folder in the repo.
These include, but are not limit to: -

- [Fragalysis Discourse Design](design_docs/Fragalysis_Discourse_v0.2.pdf)
- [Fragalysis Tags Design V1.0](design_docs/Fragalysis_Tags_Design_V1.0.pdf)
- [Fragalysis Design #651 Fix Data Download V2.0](design_docs/Fragalysis_Design_651_Fix_Data_Download_V2.0.pdf)
- [Fragalysis Job Launcher V1.0](design_docs/Fragalysis_Job_Launcher_V1.0.pdf)
- [Fragalysis Job Launcher V2.0](design_docs/Fragalysis_Job_Launcher_Phase2_V1.0.pdf)

---

[ansible]: https://github.com/ansible/ansible
[curl]: https://curl.se
[drf]: https://www.django-rest-framework.org
[environment variables]: https://docs.docker.com/compose/environment-variables/set-environment-variables/
[httpie]: https://pypi.org/project/httpie/
[jq]: https://stedolan.github.io/jq/
[kubernetes]: https://kubernetes.io
[kubernetes stack]: https://dls-fragalysis-stack-kubernetes.readthedocs.io/en/latest/index.html#
[pre-commit]: https://pre-commit.com
[readthedocs]: https://fragalysis-backend.readthedocs.io/en/latest/index.html
[sentry]: https://sentry.io/welcome/
