.. _comp-tasks:

Uploading data (Celery tasks and django template)
=================================================

Introducton
-----------
Uploading a Target or a Computed Set (see [link]) is a process that takes a long time. Typical requests handled in django
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

For example, the template for the upload page for computed sets is found at :code:`viewer/templates/viewer/upload-cset.html`.
The main parts of the template are as follows:

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

**2. Dynamic loading of upload task status and results** - this part of the template is included as a Javascript
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

The :code:`{if upload_task_id}` and :code:`{ endif }` use djangos templating language to make sure that the code
wrapped in :code:`<script>...</script>` is only executed if there is a value for :code:`upload_task_id`. The most
important variable for dynamic loading is :code:`response.data.upload_task_status` - the status of the Celery task. The
responses returned by the View are described in :ref:`Computational Data (Views) <comp-views>`

Celery task - processing and saving uploaded data
-------------------------------------------------

This task uploads a computed set from the files provided to :code:`viewer.views.UploadCSetView`, saving the data into
the relevant models specified in :ref:`Computational Data (Models) <comp-models>`. There is no separate validation
step: if anything in the files is wrong the upload is rolled back and the task fails with the error, which
:code:`viewer.views.UploadTaskView` reports back to the user.


.. autoclass:: viewer.tasks.process_compound_set
    :members:


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
