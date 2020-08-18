.. _api-intro:

RESTful API (Models, Serializers and Views)
===========================================
Introduction
------------

An API is an application program interface - a set of functions and procedures that allow the creation of applications
or scripts that can access the information, features or data of a separate operating system, application or service.

A RESTful (REpresentational State Transfer) API is an API written in a way that complies to a number of constraints to
allow it to be used as a web service. These services are written so that different services across the internet are
interoperable. RESTful web services use stateless operations -  no session information is retained by the receiver – to
allow requesting systems to manipulate text-based representations of web resources (information identifiable by URL).

RESTful APIs use HTTP methods to allow a number of different operations to be performed on a web resource over the web.
The most common operations are:

- GET: retrieve information
- POST: send information to  the server to create or update a resource
- PUT: send information to  the server to create or update a resource
- DELETE: delete a specified resource

A more complex overview of these methods (and others) is available here: https://www.w3schools.com/tags/ref_httpmethods.asp

A more thorough explaination of RESTful APIs is available here: https://searchapparchitecture.techtarget.com/definition/RESTful-API

Django REST framework (DRF)
---------------------------
In fragalysis, our RESTFul API is built using Django and Django REST framework (DRF).

DRF uses standard classes to define what information can be accessed by users, how they can access that data, and how
it is returned to the user when they make a request.

These classes first need a :code:`Model` object to operate on. :code:`Models` are the django representation of database
tables (for an example, see [link]).

Once we have a model containing data that we want to allow a user to access, we need to be able to convert that data to
a format that can be rendered in various appropriate forms (e.g. JSON or XML). To do this, we use :code:`Serializers.`
These :code:`Serializers` simply inherit a :code:`ModelSerializer` class, and all we have to tell the Serializer is what
:code:`Model` we want to operate on, and which :code:`fields` from the model we want to use (make accessible to the
user).

Once we've serialized our data, we next need to write :code:`Views`. Views define what types of requests (e.g. GET) can
be performed against a resource, and what we do with that request to return information. DRF comes with a bunch of
methods, classes and functions that allow us to do this in a simple way. We can define, using standard django methods,
what filter we want to perform on the relevant :code:`Model`, which :code:`fields` from the model we want to allow the
user to filter by (:code:`filter_fields`) and what :code:`Response` we want to post back to them.

Finally, we need to specify an :code:`endpoint` (i.e. URL) that the :code:`View` is served at, so the user can make
requests against the web service.

EXAMPLE - Model, Serializer, View and URL for Target model
----------------------------------------------------------

**Model**

The Target model contains information about a protein target. In django, we define this as an inhereted Model instance:

.. code-block:: python

    from django.db import models

    class Target(models.Model):
        """Django model to define a Target - a protein.

        Parameters
        ----------
        title: CharField
            The name of the target
        init_date: DateTimeField
            The date the target was initiated (autofield)
        project_id: ManyToManyField
            Links targets to projects for authentication
        uniprot_id: Charfield
            Optional field where a uniprot id can be stored
        metadata: FileField
            Optional file upload defining metadata about the target - can be used to add custom site labels
        zip_archive: FileField
            Link to zip file created from targets uploaded with the loader
        """
        # The title of the project_id -> userdefined
        title = models.CharField(unique=True, max_length=200)
        # The date it was made
        init_date = models.DateTimeField(auto_now_add=True)
        # A field to link projects and targets together
        project_id = models.ManyToManyField(Project)
        # Indicates the uniprot_id id for the target. Is a unique key
        uniprot_id = models.CharField(max_length=100, null=True)
        # metadatafile containing sites info for download
        metadata = models.FileField(upload_to="metadata/", null=True, max_length=255)
        # zip archive to download uploaded data from
        zip_archive = models.FileField(upload_to="archive/", null=True, max_length=255)


This model tells us which fields we want in our database table, what data type each field can contain, and some other
optional parameters, such as wether the field has to have data in it, or the maximum length of data that can be added to
the field.

**Serializer**

The Target serializer tells us what information from the Target model we want to pass back to the user when they make a
request.

.. code-block:: python

    from rest_framework import serializers
    from viewer.models import Target

    class TargetSerializer(serializers.ModelSerializer):
        template_protein = serializers.SerializerMethodField()

        def get_template_protein(self, obj):
            proteins = obj.protein_set.filter()
            for protein in proteins:
                if protein.pdb_info:
                    return protein.pdb_info.url
            return "NOT AVAILABLE"

        class Meta:
            model = Target
            fields = ("id", "title", "project_id", "protein_set", "template_protein", "metadata", "zip_archive")


The serializer uses the DRF :code:`serializers.ModelSerializer` class. We define the :code:`model` and :code:`fields` in
a :code:`Meta` subclass, where the :code:`model` is an instance of the :code:`Model` we want to operate on, and the
:code:`fields` parameter is a tuple containing the names of the fields we want to return as strings. Additionally, we
can add extra fields, and add a method to define how we get the value of the field. For example, in this
:code:`Serializer` we have added the :code:`template_protein` field, and defined how we get its value with
:code:`get_template_protein`.

**View**

This :code:`View` returns a list of information about a specific target, if you pass the :code:`title` parameter to the
request, or a list of information about all targets if you make a request against the URL.

The :code:`View` is written as a class inheriting the DRF :code:`ReadOnlyModelViewSet`, which is a standard :code:`View`
class that is read-only. That means that only GET requests can be made against this view. There are other ways to define
this for different types of view, but we won't go into detail here - this is the method we have chosen to use with most
of our standard views.

Additionally, in the actual code, you will notice that :code:`TargetView(viewsets.ReadOnlyModelViewSet)` is replaced by
:code:`TargetView(ISpyBSafeQuerySet)`. :code:`ISpyBSafeQuerySet` is a version of :code:`viewsets.ReadOnlyModelViewSet`
that includes an authentication method specific for the deployment of fragalysis at https://fragalysis.diamond.ac.uk

.. code-block:: python

    from rest_framework import viewsets
    from viewer.serializers import TargetSerializer
    from viewer.models import Target

    class TargetView(viewsets.ReadOnlyModelViewSet):
        """ DjagnoRF view to retrieve info about targets

           Methods
           -------
           url:
               api/targets
           queryset:
               `viewer.models.Target.objects.filter()`
           filter fields:
               - `viewer.models.Target.title` - ?title=<str>
           returns: JSON
               - id: id of the target object
               - title: name of the target
               - project_id: list of the ids of the projects the target is linked to
               - protein_set: list of the ids of the protein sets the target is linked to
               - template_protein: the template protein displayed in fragalysis front-end for this target
               - metadata: link to the metadata file for the target if it was uploaded
               - zip_archive: link to the zip archive of the uploaded data

           example output:

               .. code-block:: javascript

                   "results": [
                    {
                        "id": 62,
                        "title": "Mpro",
                        "project_id": [
                            2
                        ],
                        "protein_set": [
                            29281,
                            29274,
                            29259,
                            29305,
                            ...,
                        ],
                        "template_protein": "/media/pdbs/Mpro-x10417_0_apo.pdb",
                        "metadata": "http://fragalysis.diamond.ac.uk/media/metadata/metadata_2FdP5OJ.csv",
                        "zip_archive": "http://fragalysis.diamond.ac.uk/media/targets/Mpro.zip"
                    }
                ]

           """
        queryset = Target.objects.filter()
        serializer_class = TargetSerializer
        filter_permissions = "project_id"
        filter_fields = ("title",)


The docstring for this class is formatted in a way to allow a user or developer to easily read the docstring, and
understand the URL to query, how the information is queried by django, what fields can be queried against, and what
information is returned from a request against the views URL. All of the views in this documentation are written in the
same way.

**URL**

Finally, we need to define where the view is served from, in context of the root (e.g. https://fragalysis.diamond.ac.uk)
URL. The target view is served at :code:`<root>/api/targets`. In :code:`api/urls.py` we use the following lines to add
the :code:`TargetView` to that endpoint:

.. code-block:: python

    from rest_framework.routers import DefaultRouter
    from viewer import views as viewer_views

    router = DefaultRouter()
    router.register(r"targets", viewer_views.TargetView, "targets")

The DRF :code:`DefaultRouter` provides a simple, quick and consistent way of wiring ViewSet logic to a set of URLs.
Router automatically maps the incoming request to proper viewset action based on the request method type.

To make sure that we serve the URLS from :code:`api/urls.py`, we include the URLs from there in
:code:`fragalysis/urls.py`:

.. code-block:: python

    ...
    url(r"^api/", include("api.urls")),
    ...

and specify this file as the :code:`URL_ROOTCONF` in :code:`fragalysis/settings.py` - the django settings file:

.. code-block:: python

    ROOT_URLCONF = "fragalysis.urls"

If we navigate to the URL :code:`<root>/api/targets/?title=<target_name>` we are presented with the following page:

.. image:: target_api.png

This is a page automatically generated by DRF, and includes options to see what kinds of requests you can make against
this endpoint.