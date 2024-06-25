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

Basic style guidance
--------------------
All new methods (functions) in **views** _should_ inherit from classes defined in DRF.
Some methods (functions) simply accept the :code:`request`, but this is older code.

**Security**

When writing a new API endpoint, it is important to consider the security implications
of the data you are making available to the user. Much of the data in fragalysis is
_open_, i.e. available to all.

.. note::
    Open data is data associated with any :code:`Project`
    that has the :code:`open_to_public` flag set to :code:`True`.

As a general policy open data is visible to anyone, even if they are not authenticated
(logged in) but there is a policy that only logged-in users can modify or create open
data. More specifically users are required to be a member of (associated with) the
:code:`Project`` the object belongs to.

.. note::
    A user is 'associated' with a :code:`Project` (aka Proposal/Visit) if the security module in the
    project's :code:`api` package is able to find the user associated with the
    :code:`Project` by querying an external database. In diamond
    this is an **ISPyB** MySQL database external to the stack installation whose
    credentials are supplied using environment variables.

API methods that provide access to data they must ensure that the user is authenticated
and must _strive_ to ensure that the user is associated with the :code:`Project` that
the data belongs to.

In order to check whether the user has access to the object that is being created
or altered, each method must identify the :code:`Project` that the object belongs.
Given a :code:`Project` is is a relatively simple task to check that the user
has been given access to us using the security module as described above.

These actions ar simplified through the use of the :code:`ISpyBSafeQuerySet` class
to filter objects when reading and the :code:`IsObjectProposalMember` class to
check the user has access to the object when creating or altering it. These classes
rely on the definition of :code:`filter_permissions` property to direct the
search to the object's :code:`Project`.

View classes must generally inherit from :code:`ISpyBSafeQuerySet`,
which provides automatic filtering of objects. The :code:`ISpyBSafeQuerySet`
inherits from th :code:`ReadOnlyModelViewSet` view set. If a view also needs to provide
create, update or delete actions they should also inherit an appropriate
DRF **mixin**, adding support for a method so support the functionality that is
required: -

- :code:`mixins.CreateModelMixin` - when supporting objects (POST)
- :code:`mixins.UpdateModelMixin` - when supporting objects (PATCH)
- :code:`mixins.DestroyModelMixin` - when supporting delete (DELETE)

For further information refer to the `mixins`_ documentation on the DRF site.

EXAMPLE - Model, Serializer, View and URL for Target model
----------------------------------------------------------

**Model**

The Target model contains information about a protein target. In django, we define this as an inhereted Model instance:

.. code-block:: python

    from django.db import models

    class Target(models.Model):
        title = models.CharField(unique=True, max_length=200)
        init_date = models.DateTimeField(auto_now_add=True)
        project_id = models.ManyToManyField(Project)
        uniprot_id = models.CharField(max_length=100, null=True)
        metadata = models.FileField(upload_to="metadata/", null=True, max_length=255)
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

**Models**

Model definitions should avoid inline documentation, and instead use the django
:code:`help_text` parameter to provide this information. For example,
instead of doing this: -

.. code-block:: python

    class Target(models.Model):
        # The uniprot ID id for the target. A unique key
        uniprot_id = models.CharField(max_length=100, null=True)


Do this: -

.. code-block:: python

    class Target(models.Model):
        uniprot_id = models.CharField(
            max_length=100,
            null=True,
            help_text="The uniprot ID id for the target. A unique key",
        )


**View**

This :code:`View` returns a list of information about a specific target, if you pass the :code:`title` parameter to the
request, or a list of information about all targets if you make a request against the URL.

The :code:`View` is written as a class inheriting the DRF :code:`ReadOnlyModelViewSet`, which is a standard :code:`View`
class that is read-only. That means that only GET requests can be made against this view. There are other ways to define
this for different types of view, but we won't go into detail here - this is the method we have chosen to use with most
of our standard views.

Additionally, in the actual code, you will notice that :code:`TargetView(viewsets.ReadOnlyModelViewSet)` is replaced by
:code:`TargetView(ISpyBSafeQuerySet)`. :code:`ISpyBSafeQuerySet` is a version of :code:`viewsets.ReadOnlyModelViewSet`
that includes an authentication method that filters records based omn a user's
membership of the object's :code:`project`.

.. code-block:: python

    from rest_framework import viewsets
    from viewer.serializers import TargetSerializer
    from viewer.models import Target

    class TargetView(viewsets.ReadOnlyModelViewSet):
        queryset = Target.objects.filter()
        serializer_class = TargetSerializer
        filter_permissions = "project_id"
        filter_fields = ("title",)

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

.. _mixins: https://www.django-rest-framework.org/tutorial/3-class-based-views/#using-mixins
