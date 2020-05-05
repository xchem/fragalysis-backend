[![Build Status](https://travis-ci.org/xchem/fragalysis-backend.svg?branch=master)](https://travis-ci.org/xchem/fragalysis-backend)
[![stable](http://badges.github.io/stability-badges/dist/stable.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.0.1-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/xchem/fragalysis-backend/blob/master/LICENSE.txt)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis-backend.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis-backend/context:python)

# Fragalysis backend
Django server for Fragalysis with DRF API and loaders for data. Has components to serve React web-app.


# Dev environment setup
### Background

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

### Prerequisites

- Docker
- Docker-compose
- Git

## Setup

#### 1. Create project directory, e.g.

```
mkdir fragalysis
```

####2. Clone repositories inside your project's directory


You can clone original `xchem` repositories or your forked e.g. `m2ms` and checkout any branch if necessary
```
git clone https://github.com/xchem/fragalysis-backend.git
git clone https://github.com/xchem/fragalysis-frontend.git
git clone https://github.com/xchem/fragalysis-loader.git
```
Note: 
To successful build, it should exist following directories from repository cloning.
Frontend is also important, because DJANGO server will serve this directory!
```$xslt
fragalysis/fragalysis-loader/
fragalysis/fragalysis-frontend/
fragalysis/fragalysis-backend/
```

####3. Create some key data directories
In `fragalysis/` directory run script to create data directory structure
```
mkdir -p data/input/django_data/EXAMPLE
mkdir -p data/mysql/data
mkdir -p data/neo4j/data
mkdir -p data/neo4j/logs
mkdir -p data/stack/media
mkdir -p data/stack/logs
```
####4. Copy 
Copy to `fragalysis/data/input/django_data/EXAMPLE` your data, before you can launch the application.

##Start
Start `Fragalysis stack` (All infrastructure - databases + populating data)

```
docker-compose -f docker-compose.dev.yml up -d
```

`Please wait, it takes a minute until all containers are fully started.`

Test if we are running at [http://localhost:8080](http://localhost:8080)

If needed stop containers

```
docker-compose -f docker-compose.dev.yml down
```

Note: The first run will be probably not successful run `loader` and `web_dock` service, it is caused by databases
(process of first starting has long duration). Wait for successful running of `neo4j` and `mysql` services and 
retry run start script again.

##Develop mode in DJANGO
change value of `DEBUG` variable to `True` in this file
`fragalysis-backend/fragalysis/settings.py`
Note: **please don't push this change into git**