.. _comp-views:

Computational data (Views)
==========================

The models that these views work on can be found in :ref:`Computational data (Models) <comp-models>`

Where the views are used (frontend)
-----------------------------------
**Note:** For the Diamond instance of fragalysis, authentication is handled by :code:`api.security.ISpyBSafeQuerySet`,
a version of the standard DRF (see :ref:`RESTful API (Views) <api-intro>`) :code:`viewsets.ReadOnlyModelViewSet`.

- :code:`viewer.views.ComputedSetView`: This view is used to generate a list of all of the computed sets that exist for a
  specific target, so that the drop-down list for the third tab on the right-hand side (RHS) can be generated for each
  target-specific page

- :code:`viewer.views.ComputedMolAndScoreView`: This view is used to get all information about 3D molecules in a given
  computed set for a given target (computed sets obtained from :code:`viewer.views.ComputedSetView`). The smiles string
  is used to generate the 2D image shown on each molecule in the RHS computed set tab using
  :code:`viewer.views.img_from_smiles`. The scores for each molecule are displayed on the respective molecule card.

- :code:`viewer.views.cset_key`: This view is used to generate a computed set upload key at :code:`<root>/viewer/cset_key`
  , and email it to the user to allow them to upload new computed sets at :code:`<root>/viewer/upload_cset`

- :code:`viewer.views.UploadCSet`: This view is used to generate the form found at :code:`<root>/viewer/upload_cset`, and
  to pass the values to the celery tasks controlled through :code:`viewer.views.ValidateTaskView` and
  :code:`viewer.views.UploadTaskView` that validate and process the uploaded data, saving it into the relevant models
  (see :ref:`Computational data (Models) <comp-models>`)so that the data can be displayed on the RHS of a target specific
  page.

Other views listed below are currently not in use by the front-end.


View details
------------

.. autoclass:: viewer.views.ComputedSetView
    :members:

.. autoclass:: viewer.views.ComputedMoleculesView
    :members:

.. autoclass:: viewer.views.NumericalScoresView
    :members:

.. autoclass:: viewer.views.TextScoresView
    :members:

.. autoclass:: viewer.views.CompoundScoresView
    :members:

.. autoclass:: viewer.views.ComputedMolAndScoreView
    :members:

.. autoclass:: viewer.views.cset_key
    :members:

.. autoclass:: viewer.views.UploadCSet
    :members:

.. autoclass:: viewer.views.UploadTSet
    :members:

.. autoclass:: viewer.views.ValidateTaskView
    :members:

.. autoclass:: viewer.views.UploadTaskView
    :members:

.. autoclass:: viewer.views.cset_download
    :members:

.. autoclass:: viewer.views.pset_download
    :members:

