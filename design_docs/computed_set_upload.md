# Compound Set Upload Specification

Computed sets of compounds can be loaded into Fragalysis (aka the right hand side).

The original specification for this can be found
[here](https://discuss.postera.ai/t/providing-computed-poses-for-others-to-look-at/1155/8).
This document aims to provide a succinct summary of the current status in Fragalysis.

*Note that the current upload functionality is currently in flux so these details will probably change.*

## Current Specification - ver_1.2

The upload format for compounds will be allowed in one of two ways:

1. A single sdf file
2. A single sdf file plus pdb files for the ligands to be loaded into in fragalysis.

The sdf files for these two options will have a standardised format, to allow the following options:

- The fragments that inspired the design of each molecule can be specified
- The protein (in pdb file format) for each molecule can be specified
- Any number of ‘properties’ or ‘scores’ can be specified.

### The SD-file format is as follows:

The sdf file name will be: compound-set_<name>.sdf with <name> replaced with the name you wish to give it. e.g. compound-set_fragmenstein.sdf

A **blank** molecule will be the first in the sdf:

- This molecule will contain all of the same fields as the sdf, containing a description of those fields.
- The 3D coordinates of this molecule can be anything - they will be ignored.
- The name (title line) of this molecule should be the file format specification version e.g. ver_1.2 (as defined in this document)
- The molecule should have the following compulsory fields:
- ref_url - the url to the forum post that describes the work
- submitter_name - the name of the person submitting the compounds
- submitter_email - the email address of the submitter
- submitter_institution - the submitters institution
- generation_date - the date that the file was generated in format yyyy-mm-dd
- method - a name for the method used to generate the compound poses

**NB: all of the compulsory fields for the blank mol can be included for the other molecules, but they will be ignored**

**Every other molecule** in the sdf file will be assumed to be a molecule that is a computed molecule, and should:

- Have the same properties as the blank molecule, but with their values instead of description. - for the ref_url field, you can leave this blank, as it will be ignored for molecules that are not the blank molecule
- Have a name that is meaningful, and will eventaully be displayed in Fragalysis - Use the PostEra submission ID, or if not available, a name that is meaningful (e.g. the PDB code, name used in publication, etc.)
- Have the three following compulsary property fields:
- ref_mols - a comma separated list of the fragments that inspired the design of the new molecule (codes as they appear in fragalysis - e.g. x0104_0,x0692_0)
- ref_pdb - either:
    - the file path of the pdb file in the uploaded zip file:
      Example: If you upload a file called references.zip that contains a pdb file new_protein.pdb, the corresponding path in the ref_pdb file would be references/new_protein.pdb
    - the code to the fragment pdb from Fragalysis that should be used (e.g. x0692_0)
- original SMILES - the original smiles of the compound before any computation was carried out

**NB: only properties with numerical (or boolean) values will be displayed in Fragalysis - this will be reviewed at a later date**
