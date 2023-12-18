.. _stack-setup:

Local developer environment setup
=================================
Background
----------

The stack consists of three services, running as containers: -

- a MySQL database
- a neo4j graph database
- the fraglaysis stack
- a transient data loader container

The stack is formed from code resident in a number of repositories.
Begin by forking repositories you anticipate editing (although you really want
to consider forking all the repositories as this is a relatively low-cost
operation).

The repositories are:

- [xchem/fragalysis](https://github.com/xchem/fragalysis)
- [xchem/fragalysis-frontend](https://github.com/xchem/fragalysis-frontend)
- [xchem/fragalysis-backend](https://github.com/xchem/fragalysis-backend)
- [xchem/fragalysis-stack](https://github.com/xchem/fragalysis-stack)
- [xchem/fragalysis-loader](https://github.com/xchem/fragalysis-loader)

Prerequisites
-------------

- Docker
- Docker-compose
- Git
- NodeJS (v12)
- Yarn
- Some target data

Setup
-----

**1. Create project directory, e.g.**

.. code-block:: bash

    mkdir fragalysis

**2. Clone repositories inside your project's directory**


You can clone original :code:`xchem` repositories or your forked e.g. :code:`m2ms` and checkout any branch if necessary

.. code-block:: bash

    git clone https://github.com/xchem/fragalysis-backend.git
    git clone https://github.com/xchem/fragalysis-frontend.git
    git clone https://github.com/xchem/fragalysis-loader.git
    git clone https://github.com/InformaticsMatters/dls-fragalysis-stack-openshift.git


**3. Create some key data directories**


In the :code:`fragalysis/` directory run script to create data directory structure

.. code-block:: bash

    mkdir -p data/input/django_data/EXAMPLE
    mkdir -p data/neo4j/data
    mkdir -p data/neo4j/logs
    mkdir -p data/stack/media
    mkdir -p data/stack/logs
    mkdir -p data/media/computed_set_data
    mkdir -p data/postgre/data

**4. Build images locally**

*Optional*


.. code-block:: bash

    pushd fragalysis-backend || exit
    docker build . -t xchem/fragalysis-backend:latest
    popd || exit

    pushd fragalysis-loader || exit
    docker build . -t xchem/fragalysis-loader:latest
    popd || exit

    pushd fragalysis-stack || exit
    docker build . -t xchem/fragalysis-stack:latest
    popd || exit

*Mandatory*

.. code-block:: bash

    pushd dls-fragalysis-stack-openshift/images/loader || exit
    docker build . -f Dockerfile-local -t loader:latest
    popd || exit

*Optional*

.. code-block:: bash

    pushd dls-fragalysis-stack-openshift/images/graph || exit
    docker build . -t xchem/graph:latest
    popd || exit


**4. Populating the database**

Copy your example data to :code:`fragalysis/data/input/django_data/EXAMPLE`, before you can launch the application.

Launch the stack
-----------------

1. To build the :code:`Fragalysis stack` (All infrastructure - databases + populating data):

.. code-block:: bash

    cd fragalysis-backend
    docker-compose -f docker-compose.dev.yml up -d

*Note*: sometimes the database does not finish building before the data loader starts. If this happens, allow the build
to finish, kill it, and then run the build again

2. Navigate to :code:`localhost:8080` to check that the front-end is running, and that the data has successfully loaded
   into the back-end


3. If needed, stop containers

.. code-block:: bash

    docker-compose -f docker-compose.dev.yml down

**Note:** If you want to connect to get into the stack container (:code:`web_dock`) run:

.. code-block:: bash

    docker exec -it web_dock /bin/bash

Troubleshooting
---------------
**Problem:** The database keeps throwing errors when I try to build the stack

**Problem:** The data in my local stack has disappeared.

**Solution:** Sometimes the postgres database becomes corrupted, or there are problems with migrations. To fix:

- Delete the 'mysql_data' folder located in the data folder (Step 3)
- Delete the migrations by running the commands below in the project directory (Step 1)

.. code-block:: bash

    find . -path "*/migrations/*.py" -not -name "__init__.py" -delete
    find . -path "*/migrations/*.pyc"  -delete

- In the 'fragalysis-frontend' folder rerun the docker compose command another two times

.. code-block:: bash

    docker-compose -f docker-compose.dev.yml up

----

**Problem:** I keep getting 500 codes when I make changes to the backend code, but I don't know why

**Solution:** Turn on Django's debug mode

- change the value of :code:`DEBUG` to :code:`True` in :code:`fragalysis-backend/fragalysis/settings.py`

Note: **please don't push this change into git**

----

**Problem:** No matter what I try, the stack won't build anymore (it did before)

**Solution:** Clean up and rebuild:

.. code-block:: bash

    # remove all old images
    docker system prune -a

    # rebuild the loader image (and any other images needed)
    pushd dls-fragalysis-stack-openshift/images/loader || exit
    docker build . -f Dockerfile-local -t loader:latest
    popd || exit #etc

    # re-run docker compose
    docker-compose -f docker-compose.dev.yml up
