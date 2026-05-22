Miscellaneous data
==================

Miscellaneous Views
-------------------

.. autofunction:: viewer.views.img_from_smiles

.. autofunction:: viewer.views.highlight_mol_diff

.. autofunction:: viewer.views.get_open_targets

Tags
====

Tags provide rich, structured labelling of data, with the following functionality:

- Tags can be attached to Site Observations and Session Projects. The functionality uses a base :code:`Tag` model class with sub-classes (:code:`SiteObservationTag` and :code:`SessionProjectTag`) containing the fields relevant to each.

- Tag information includes a Discourse URL and a JSON field containing any additional information.

- Tags are attached to Categories so they can easily be distinguished on the Front End. The category table is loaded with categories such as "Sites", "Series", "Forum" and "Other".

- There are supporting APIs to allow access to the :code:`SessionProjectTag` and :code:`SiteObservationTag` models.



Tag Model details
-----------------

.. autoclass:: viewer.models.TagCategory

.. autoclass:: viewer.models.Tag

.. autoclass:: viewer.models.SiteObservationTag

.. autoclass:: viewer.models.SessionProjectTag

Tag Views
---------


.. autoclass:: viewer.views.TagCategoryView
    :members:


.. autoclass:: viewer.views.SiteObservationTagView
    :members:


.. autoclass:: viewer.views.SessionProjectTagView
    :members:
