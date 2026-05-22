Serving static files (media)
============================

There are some files, such as PDB files, that need to be served by fragalysis. These are either:

- **Served to the front-end** - files retrieved from a specific URL to be loaded and rendered in the front-end - most
  commonly PDB files loaded by NGL into the central 3D viewer

- **Served as downloads** - files retrieved from a specific URL to be downloaded directly by a user.

It is not recommended to use Django to serve static files in production, so fragalysis makes use of nginx to serve its
media.
(https://docs.djangoproject.com/en/3.1/howto/static-files/deployment/#serving-static-files-from-a-dedicated-server).

NGINX
-----

According to the NGINX website (https://www.nginx.com/resources/glossary/nginx/):

    NGINX is open source software for web serving, reverse proxying, caching, load balancing, media streaming, and more.
    It started out as a web server designed for maximum performance and stability.

Adding an endpoint that serves media
------------------------------------

In fragalysis, we use a :code:`FileField` in a :code:`Model` to specify a file, and where it is stored. For example, in
the target model we have two files that are served as files:

.. code-block:: python

    class Target(models.Model):
        title = models.CharField(max_length=200, help_text="A title, i.e. Mpro")
        init_date = models.DateTimeField(auto_now_add=True)
        project = models.ForeignKey(Project, on_delete=models.CASCADE)
        uniprot_id = models.CharField(
            max_length=100,
            null=True,
            help_text="The uniprot ID id for the target. A unique key",
        )
        metadata = models.FileField(
            upload_to="metadata/",
            null=True,
            max_length=255,
            help_text="Optional file upload defining metadata about the target."
            " Can be used to add custom site labels",
        )
        zip_archive = models.FileField(
            upload_to="archive/",
            null=True,
            max_length=255,
            help_text="Link to zip file created from targets uploaded with the loader",
        )

The :code:`metadata` file and the :code:`zip_archive` file. To make sure that these :code:`Files` are served as download links
, we have to tell nginx to allow access to the media area where the files are uploaded/saved to.

We define this in :code:`django_nginx.conf`, for example for the :code:`metadata` file:

.. code-block:: bash

    location /metadata/ {
            alias /code/media/metadata/;
            internal;
        }

This tells nginx that the url :code:`/metadata/` can be used to access files hosted on the back-end under
:code:`/code/media/metadata`.

In order to tell django to serve these files at the given URL, we need a :code:`View`, just like any other endpoint. The
view for serving :code:`metadata` files can be found in :code:`media_serve/views.py`:

.. code-block:: python

    from api.security import ISPyBSafeStaticFiles
    from viewer.models import Target

    def metadata_download(request, file_path):
        """
        Download a metadata file by nginx redirect
        :param request: the initial request
        :param file_path: the file path we're getting from the static
        :return: the response (a redirect to nginx internal)
        """
        ispy_b_static = ISPyBSafeStaticFiles()
        ispy_b_static.model = Target
        ispy_b_static.request = request
        ispy_b_static.permission_string = "project"
        ispy_b_static.field_name = "metadata"
        ispy_b_static.content_type = "application/x-pilot"
        ispy_b_static.prefix = "/metadata/"
        ispy_b_static.input_string = file_path
        return ispy_b_static.get_response()


For our files, we're authenticating user access to file downloads with :code:`ISPyBSafeStaticFiles`, which is a custom
class inheriting :code:`ISPyBSafeQuerySet(viewsets.ReadOnlyModelViewSet)`: the same View method we use to authenticate
user access for all non-media views. The difference between :code:`ISPyBSafeQuerySet` and :code:`ISPyBSafeStaticFiles`
is that :code:`ISPyBSafeStaticFiles` contains a method that sets the context of the response using NGINX's redirect
method, returning the response including the file as an attachment:

.. code-block:: python

    class ISPyBSafeStaticFiles:

        def get_queryset(self):
            query = ISPyBSafeQuerySet()
            query.request = self.request
            query.filter_permissions = self.permission_string
            query.queryset = self.model.objects.filter()
            queryset = query.get_queryset()
            return queryset

        def get_response(self):
            try:
                queryset = self.get_queryset()
                filter_dict = {self.field_name + "__endswith": self.input_string}
                object = queryset.get(**filter_dict)
                file_name = os.path.basename(str(getattr(object, self.field_name)))

                if hasattr(self, 'file_format'):
                    if self.file_format=='raw':
                        file_field = getattr(object, self.field_name)
                        filepath = file_field.path
                        zip_file = open(filepath, 'rb')
                        response = HttpResponse(FileWrapper(zip_file), content_type='application/zip')
                        response['Content-Disposition'] = 'attachment; filename="%s"' % file_name

                else:
                    response = HttpResponse()
                    response["Content-Type"] = self.content_type
                    response["X-Accel-Redirect"] = self.prefix + file_name
                    response["Content-Disposition"] = "attachment;filename=" + file_name

                return response
            except Exception:
                raise Http404

Finally, just as with any other view, we have to specify a url. For example, for the :code:`metadata` file, we specify
in :code:`media_serve/urls.py`:

.. code-block:: python

    path("metadata/<file_path>/", views.metadata_download, name="get_metadata"),
