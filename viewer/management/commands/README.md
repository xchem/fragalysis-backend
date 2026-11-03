# Commands for managing specific tasks in Fragalysis

The commands in this directory can be run from the container directly using
a command of the form:
```
python manage.py <python script> <options>
```

It should be possible to run them either manually or as part of a cron job to
do maintenance and data intergrity-style tasks.

It would be good if all commands are lightly documented here.

## start_service_queries

Defined in app `service_status`.

```
python manage.py start_service_queries
```

Starts the external service health check tasks. Reads the required
tasks from `ENABLE_SERVICE_STATUS` environment variable and launches
the corresponding task (test function must be defined in
`service_status/services.py`).


## services

Defined in app `service_status`.

```
python manage.py services --enable [comma separated service names]

# or

python manage.py services --disable [comma separated service names]
```

Starts or stops service health check queries individually.
