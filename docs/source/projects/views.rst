.. _proj-views:

Project data (Views)
==========================

The models that these views work on can be found in :ref:`Project data (Models) <proj-models>`

Where the views are used (frontend)
-----------------------------------

- :code:`viewer.views.SnapshotsView`: the json data in the :code:`data` field in this view is used by the reducers on
  the front-end to re-create the context of any page that a user saves to a project, or saves as a standalone snapshot.
  The PUT method of this view is what is used by the front-end button push of 'save' or 'share' in the menu bar of a
  target specific page to send the state of the page to the view so that the data is saved and accessed later. The data
  can be accessed either through an endpoint, or if the user is logged in when they save the session, as part of its
  associated project, all of which are listed on the landing page when a user is logged in. If the user is logged in and
  navigates to a specific project, the project navigator will appear on the bottom right-hand side as a tree, which the
  user can navigate through to access different project snapshots. Information such as who created the snapshot, its
  title and description are also displayed here.

- :code:`viewer.views.SessionProjectsView`: This view is used to create a project (collection of snaphots) when a user
  is logged in and clicks the '+' button on the projects pane of the landing page. The GET method of this view is used
  to retrieve information about user-created projects, and display it in the projects pane of the landing page.

There are also two subsidiary views for Snapshot and Session Project to capture and store actions by the user on the
frontend. These allow the user to recreate the current state and allow the user to go backwards and forwards through
previous actions. These views are:

- :code:`viewer.views.SnapshotActionsView`
- :code:`viewer.views.SessionActionsView`

Discourse
---------

If a Discourse platform has been configured and a Discourse API key is supplied in the environmental parameters,
then the :code:`DiscoursePostView` can be used to create sub-categories, topics and posts. The intention is that:
- Subcategories will be created when a Target is uploaded (with the name of the Target)
- Creating a Session Project for the Target will also allow a new topic within the sub-category (Target)
- Creating a Snapshot within a Session Project will also allow a new post to be created within the topic

In all cases a URL is returned allowing the user to immediately navigate to the post if desired.


View details
------------

.. autoclass:: viewer.views.SnapshotsView
    :members:

.. autoclass:: viewer.views.SnapshotActionsView
    :members:

.. autoclass:: viewer.views.SessionProjectsView
    :members:

.. autoclass:: viewer.views.SessionActionsView
    :members:

.. autoclass:: viewer.views.DiscoursePostView
    :members:
