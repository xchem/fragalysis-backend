Download data
=============

Fragalysis provides functionality to flexibly download subsets of data as follows:

- **Subset of Target Data** - Given a permitted target, a list of protein codes, and a series of booleans indicating which
  information is requested, the DownloadStructures view will construct a zip file of the requested information (if available).
  For example, the user could request to download only PDB files for a particular list of protein codes.

- **Subset of Computed Set Data** - Constructs a csv file for download based on a dictionary constructed in the react
  front end from the computed sets.


 Views
---------


.. autoclass:: viewer.views.DownloadStructures
    :members:

.. autoclass:: viewer.views.DictToCsv
    :members:
