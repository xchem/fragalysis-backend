Miscellaneous data
==================

Miscellaneous Views
-------------------

.. autoclass:: viewer.views.img_from_smiles
    :members:

.. autoclass:: viewer.views.highlight_mol_diff
    :members:

.. autoclass:: viewer.views.get_open_targets
    :members:

Tags
====

From April 2021, the Tag functionality is being upgraded. Before, tags were simple strings added to a
Session Project. Following this introduction, the tag information will be much richer, providing (eventually) the
following new functionality:

- Tags can be attached to Molecules and Session Projects. The functionality uses a base "Tag" model class with sub-classes containing the fields relevant to Molecules and Session Projects.

- Tag information has been expanded to include a Discourse URL and a JSON field containing any additional information.

- Tags are now attached to Categories so they can easily be distinguished on the Front End. Although this is expandable the category table will be initially loaded with "Sites", "Series", Forum" and "Other".

- The tags relating to "Sites" will replace the current Site functionality. In the cross over phase, both the old sites (Molgroups) and new Site Tags will be active and loaded each time a new target is uploaded. The intention is to expand the Site Tag functionality in the future when the front end has been changed and the Molgroup functionality is no longer used.

- A new API has been introduced to load all Target and Molecule/Tag information for the initial load of a Target. This is intended to replace the current functionality that loads a dataset molecule by molecule.

- There are supporting APIs to allow access to the new SessionProjectTag, MoleculeTag models.



Tag Model details
-----------------

.. autoclass:: viewer.models.TagCategory

.. autoclass:: viewer.models.Tag

.. autoclass:: viewer.models.MoleculeTag

.. autoclass:: viewer.models.SessionProjectTag

Tag Views
---------


.. autoclass:: viewer.views.TagCategoryView
    :members:


.. autoclass:: viewer.views.MoleculeTagView
    :members:


.. autoclass:: viewer.views.SessionProjectTagView
    :members:


.. autoclass:: viewer.views.TargetMoleculesView
    :members:
