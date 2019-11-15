# Fragalysis REST Endpoints

The Fragalysis backend is Python Django application serving up data using a REST API.
The [django-rest-swagger](https://django-rest-swagger.readthedocs.io/en/latest/) module
is enabled which serves up Swagger (OpenAPI) descriptions of the endpoints.

These are some key endpoints for Fragalysis.
They are defined in the various views.py files, with the key one being
[fragalysis/views.py]().

| Relative URL  | Prod site URL                                  | Description |
| ------------- | ---------------------------------------------- | ----------- |
| /api/         | https://fragalysis.diamond.ac.uk/api/          | Django API endpoints |
| /api/swagger/ | https://fragalysis.diamond.ac.uk/api/swagger/  | API endpoints as Swagger | 

**NOTE:** the Swagger implementation is currently not working fully as some endpoints are
using HTTP not HTTPS which blocks the Swagger web UI. 

Examples:
* Get all targets: https://fragalysis.diamond.ac.uk/api/targets/
* Get a particular target: https://fragalysis.diamond.ac.uk/api/targets/2/

## TODO

* Get Swagger working fully
* Fully document API endpoints and their parameters 
* Describe how to authenticate
