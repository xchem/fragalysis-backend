# Download Structures Notes

Firstly, note that everything comes from the aligned path.
The original plan was for the contents of the zip file to come from the links that are currently stored in the Fragalysis database, so I did the design based on this.

The extra work for the crystalographic path is detailed in https://github.com/m2ms/fragalysis-frontend/issues/673. This needs updating I think based on what we know now.

The mapping for the aligned files stored in the database (the ones we currently output) is as follows:

| boolean in api      | File template for source in aligned directory                                                      |
|---------------------|----------------------------------------------------------------------------------------------------|
| 'pdb_info'          | "_apo.pdb"                                                                                         |
| 'bound_info'        | "_bound.pdb"                                                                                       |
| 'cif_info'          | There is currently no ‘CIF’ mapping so this will be empty for current uploads.                     |
| 'mtz_info'          | ".mtz"                                                                                             |
| 'map_info'          | The ‘PMAP’ mapping (was "_event.ccp4") is commented out so this will be empty for current uploads. |
| 'sigmaa_info'       | "_2fofc.map"                                                                                       |
| 'diff_info'         | "_fofc.map"                                                                                        |
| 'event_info'        | "_event.ccp4"                                                                                      |
| 'trans_matrix_info' | "_transform.json"                                                                                  |
| ‘sdf_info’          | ".sdf"                                                                                             |

The boolean names are identical to the field names in the viewer_protein and viewer_molecule tables.
