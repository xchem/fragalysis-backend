.. _comp-tasks:

Uploading computational data (Celery tasks and django template)
===============================================================

Introducton
-----------
Uploading a Computed Set (see [link]) is a process that takes a long time. Typical requests handled in django
applications usually only take milliseconds. Because this process takes a long time, we have implemented the required
components as Celery tasks, and written views that can handle kicking off the relevant tasks, and retrieve the results
from the tasks by pinging the message queue that brokers the tasks (Redis).

Celery is a task queuing software package, allowing execution of asynchronous workloads. A synchronous operation blocks
a process till the operation completes. An asynchronous operation is non-blocking and only initiates the operation. In
the context of a web-server, this means we can kick off an operation through an endpoint, and query another endpoint to
periodically check if the task has completed, and when it has, retrieve the results.

This is advantageous when using processes that can take a long time, because browsers often have built in settings where
the longest it will wait for the result of a request is around 30 seconds. If the task we want to run takes longer than
that, the user will receive a timeout message.

Instead, we can use a django template that uses a bit of javascript to dynamically ping the endpoint that checks on the
status of a job, and renders a result or message in the same page when it is completed or fails.

A good tutorial about using Celery and Redis with django can be found here:
https://stackabuse.com/asynchronous-tasks-in-django-with-redis-and-celery/


Template - the upload page for computed sets
--------------------------------------------
A django template (https://docs.djangoproject.com/en/3.1/topics/templates/) is an HTML document that can dynamically
load content from django objects, such as models and views. For example, we could query the Target model, and dynamically
create a list of all Target names in a template.

The template for the upload page for computed sets is found at :code:`viewer/templates/viewer/upload-cset.html`. The
main parts of the template are as follows:

**1. The upload form** - this part of the template makes use of :code:`viewer.forms.CSetForm`: a version of a :code:`Model`
that is used to describe what information can be posted as a request through a form contained in a template.

.. autoclass:: viewer.forms.CSetForm
    :members:

The code that shows this form is:

.. code-block:: html


    <form method="post" enctype="multipart/form-data">
        {% csrf_token %}
        {{ form.as_ul }}
        <button type="submit">Submit</button>
    </form>

This code uses standard html tags to show that we're using a form, and uses django's templating language to insert the
form as a list (:code:`{{ form.as_ul}}`)

**2. Dynamic loading of validate task status and results** - this part of the template is included as a Javascript
script, and uses a function to periodically ping the task endpoint, updating what is displayed on the page depending on
the tasks status:

.. code-block:: html

    {% if validate_task_id %}
      <script>
      var content = document.getElementById('content');
      content.innerHTML = "";
      var taskUrl = "{% url 'validate_task' validate_task_id=validate_task_id %}";
      var dots = 1;
      var progressTitle = document.getElementById('progress-title');
      updateProgressTitle();
      var timer = setInterval(function() {
        updateProgressTitle();
        axios.get(taskUrl)
          .then(function(response){
            var taskStatus = response.data.validate_task_status
            if (taskStatus === 'SUCCESS') {

              var content = document.getElementById('links');
              content.innerHTML = response.data.html;

              clearTimer('');
            }

            else if (taskStatus === 'FAILURE') {

                clearTimer('An error occurred - see traceback below');
                var content = document.getElementById('links');
                content.innerHTML = response.data.validate_traceback;

            }
          })
      }, 800);

      function updateProgressTitle() {
        dots++;
        if (dots > 3) {
          dots = 1;
        }
        progressTitle.innerHTML = 'validating files';
        for (var i = 0; i < dots; i++) {
          progressTitle.innerHTML += '.';
        }
      }
      function clearTimer(message) {
        clearInterval(timer);
        progressTitle.innerHTML = message;
      }
     </script>
    {% endif %}

The :code:`{if validate_task_id}` and :code:`{ endif }` use djangos templating language to make sure that the code
wrapped in :code:`<script>...</script>` is only executed if there is a value for :code:`validate_task_id`. This value
comes from the Validate task, which is described below.

The different variables (:code:`var:...`) are used to decide what to render on the page in the different elements
defined in functions (e.g. :code:`document.getElementById('links')` - which looks for the HTML div named :code:`links`).

The most important variable for dynamic loading is the :code:`response.data.validate_task_status` variable - this is the
status of the validate task from the Celery task, returned by the Validate task, which is described below.

The responses returned by the View are described in :ref:`Computational Data (Views) <comp-views>`

**3. Dynamic loading of upload task status and results** - this part of the template is included as a Javascript
script, and uses a function to periodically ping the task endpoint, updating what is displayed on the page depending on
the tasks status:

.. code-block:: html

      {% if upload_task_id %}
          <script>
              var content = document.getElementById('content');
              content.innerHTML = "";
              var taskUrl = "{% url 'upload_task' upload_task_id=upload_task_id %}";
              var dots = 1;
              var progressTitle = document.getElementById('progress-title');
              updateProgressTitle();
              var timer = setInterval(function() {
                updateProgressTitle();
                axios.get(taskUrl)
                  .then(function(response) {
                      var taskStatus = response.data.upload_task_status
                      if (taskStatus === 'SUCCESS') {
                          var validatedStatus = response.data.validated
                          if (validatedStatus === 'Not validated') {

                              var content = document.getElementById('links');
                              content.innerHTML = response.data.html;

                              clearTimer('');

                          }
                          if (validatedStatus === 'Validated') {
                              clearTimer('Your files were uploaded! The download links are:');

                              var url_a = response.data.results.cset_download_url;
                              var content = document.getElementById('links');
                              var a = document.createElement("a");
                              var link = document.createTextNode("    Compound Set    ");
                              a.appendChild(link);
                              a.title = 'compound set';
                              a.href = url_a;
                              content.appendChild(a);

                              var br = document.createElement('br');
                              content.appendChild(br);

                              var url_b = response.data.results.pset_download_url;
                              var b = document.createElement("a");
                              var link_b = document.createTextNode("    Protein Set    ");
                              b.appendChild(link_b);
                              b.title = 'protein set';
                              b.href = url_b;
                              content.appendChild(b);

                          }
                          var moleculesProcessed = response.data.processed

                          if (moleculesProcessed === 'None') {

                              var content = document.getElementById('links');
                              content.innerHTML = response.data.html;

                              clearTimer('');

                          }
                      }
                        else if (taskStatus === 'FAILURE') {

                            clearTimer('An error occurred - see traceback below');
                            var content = document.getElementById('links');
                            content.innerHTML = response.data.upload_traceback;

                    }
                  })
              }, 800);

              function updateProgressTitle() {
                dots++;
                if (dots > 3) {
                  dots = 1;
                }
                progressTitle.innerHTML = 'processing uploaded files';
                for (var i = 0; i < dots; i++) {
                  progressTitle.innerHTML += '.';
                }
              }
              function clearTimer(message) {
                clearInterval(timer);
                progressTitle.innerHTML = message;
              }
         </script>
      {% endif %}

This code works in the same way as the Javascript code for the validate task, but instead uses the Upload task described
in :ref:`Computational Data (Views) <comp-views>`

Celery task - validating uploaded data
--------------------------------------

The first task that has to be completed when uploading a computed set is validation of the data. This task checks the
format of the uploaded SDF file provided to :code:`viewer.views.UploadCSetView` to make sure it is in the correct format
and contains all of the required information (specified here: [link]) to upload save the data to the database through
the :code:`viewer.views.UploadTaskView` into the relevant models specified in
:ref:`Computational Data (Models) <comp-models>`.


.. autoclass:: viewer.tasks.validate_compound_set
    :members:


Celery task - processing and saving uploaded data
-------------------------------------------------

The second task that has to be completed when uploading a computed set is the upload itself. This task checks takes the
output of :code:`viewer.tasks.validate` - the uploaded files must be validated before their data can be saved to the
database.


.. autoclass:: viewer.tasks.process_compound_set
    :members:


Uploading Target data sets
--------------------------

From 2021, a similar process has been developed for Target data sets to replace the existing process
that are created from the source data using the Fragalysis Loader repo.
Target data sets are initially created in the same way as currently using the Fragalysis api repo (
found here: https://github.com/xchem/fragalysis-api ).

The template for the upload page for target sets is found at :code:`viewer/templates/viewer/upload-tset.html`. The
main parts of the template are similar to the computed sets.

The view controlling interaction is: :code:`viewer.views.UploadTSetView`.

The upload form is:

.. autoclass:: viewer.forms.TSetForm
    :members:

Note that the user must be **logged on to Fragalysis** to be able to use the form. When the user is logged on, the
email field is filled with the contents of the user.email field. This is used to send a notification to the user
on completion of an validation/upload task (see below).

Dynamic loading of validate task status and results otherwise works in a similar way to the computational sets.

Celery task - validating Target data
--------------------------------------

The first task that has to be completed when uploading a target set is validation of the data. This task checks the
format of the folder structure provided to :code:`viewer.views.UploadTSetView` to make sure it is in the correct format.

.. autoclass:: viewer.tasks.validate_target_set
    :members:


Celery task - processing Target data
-------------------------------------------------

The second task that has to be completed when uploading a computed set is the upload itself. This task checks takes the
output of :code:`viewer.tasks.validate_target_set` - the uploaded files must be validated before their data can be
saved to the database through the :code:`viewer.views.UploadTaskView` into the relevant models specified in
:ref:`Crystallographic data (Models) <crys-models>`.

.. autoclass:: viewer.tasks.process_target_set
    :members:

