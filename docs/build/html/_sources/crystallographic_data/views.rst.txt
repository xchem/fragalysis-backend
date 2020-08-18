.. _crys-views:

Crystallographic data (Views)
==========================

The models that these views work on can be found in :ref:`Crystallographic data (Models) <crys-models>`

Where the views are used (frontend)
-----------------------------------
**Note:** For the Diamond instance of fragalysis, authentication is handled by :code:`api.security.ISpyBSafeQuerySet`,
a version of the standard DRF (see :ref:`RESTful API (Views) <api-intro>`) :code:`viewsets.ReadOnlyModelViewSet`.

- :code:`viewer.views.TargetView`: This view is used on the landing page
  (e.g. https://fragalysis.diamond.ac.uk/viewer/react/landing) to populate the Target list with either open targets (if
  the user is not logged in) or the Target list of targets that the user is authenticated to see.

- :code:`viewer.views.ProteinView`: This view is used on each target specific page to populate protein information on the
  left-hand side (LHS)

- :code:`viewer.views.ProteinPDBInfoView`: This view is used to retrieve the raw pdb file info to display in the central
  3D viewer when a protein for a specific ligand is turned on from the LHS when the 'P' button is clicked in the hit
  navigator.

- :code:`viewer.views.ProteinPDBBoundInfoView`: This view was previously used to store the bound-state information for
  each crystal in a target set, and the full list for each target was zipped into a file by the front-end for download
  with the download button in the menu bar of a target specific page. This has been replaced by serving a zip file of all
  of the information uploaded to fragalysis for a specific target. (see: [link])

- :code:`viewer.views.ProteinMapInfoView`: This view is currently un-used, but is intended to serve the raw electron
  density data for each ligand when the 'D' button is clicked for a ligand in the hit navigator

- :code:`viewer.views.MoleculeView`: This view is used to retrieve 2D information about molecules - the properties in the
  view are displayed on each ligand card on the LHS.

- :code:`viewer.views.MolImageView`: This view is used to retreive the 2D image created by the loader, and display it on
  the LHS

- :code:`viewer.views.CompoundView`: This view is used to obtain the 3D coordinates for a molecule, and display them in
  the central 3D viewer on a target specific page when the 'L' button on a specific ligand card in the hit navigator is
  turned on.


Other views listed below are currently not in use by the front-end.


View details
------------

.. autoclass:: viewer.views.TargetView
    :members:

.. autoclass:: viewer.views.ProteinView
    :members:

.. autoclass:: viewer.views.ProteinPDBInfoView
    :members:

.. autoclass:: viewer.views.ProteinPDBBoundInfoView
    :members:

.. autoclass:: viewer.views.ProteinMapInfoView
    :members:

.. autoclass:: viewer.views.MoleculeView
    :members:

.. autoclass:: viewer.views.MolImageView
    :members:

.. autoclass:: viewer.views.CompoundView
    :members:

.. autoclass:: viewer.views.CompoundImageView
    :members:

