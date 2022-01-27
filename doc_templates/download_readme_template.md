# Fragalysis Download Directory Sturcture

A fragalysis download will contain 2-3 folders and some additional files at the top level directory.

At the top level there are 2 files. `metadata.csv` and `smiles.smi`. These are both plain-text files. `metadata.csv` will contain information about the context of each ligand and may provide a convenient way to browse through smiles, site labels and PDB codes for each ligand. `smiles.smi` contains a list of all smiles strings that you have downloaded separated by commans.

`[TARGETNAME]_combined.sdf` may also be present which will contain all the ligand sdf files in a single sdf file.

## aligned directory

The aligned directory contains a subdirectory for each ligand that was selected for downloading.

### Contents of aligned ligand subdirectory

Depending on your selection of options when downloading the data the follow file suffixes may be present

- [ligand_name]\_apo.pdb --- protein model without ligand bound
- [ligand_name]\_bound.pdb - protein model with ligand bound
- [ligand_name]\_event.(map/ccp4) - Event Electron density cut to around 12 Angstrom around the ligand. This has a higher signal-to-noise ratio which will amplify the evidence of ligand occupancy
- [ligand_name]\_2fofc.(map/ccp4) - estimate of the true electron density from diffraction data and atomic model. Cut to around 12 Angstrom around the ligand.
- [ligand_name]\_fofc.(map/ccp4) - difference electron density map, negative density typically represents where no electron density is found but exists in the atom model. Positive densities represent electron density without mapped atom model. Cut to around 12 Angstrom around the ligand.
- [ligand_name].sdf - The Ligand molecule in sdf format
- [ligand_name]\_transform.json - Tranformation matrix and vector in json format used to align all data together.

## crystallographic directory

The crystallographic folder contains the unprocessed versions of all data found in the aligned folder. As one crystal can have mutliple ligands we provide the input crystallographic files once to avoid redundancy and keep download sizes to a minimum.

### Contents of crystal subdirectory

Depending on your selection of options when downloading the data the follow file suffixes may be present:

- [crystal_name].pdb
- [crystal_name].mtz Reflection data corresponding to pdb file.
- [crystal_name]\_event.mtz Event Backgroud corrected reflection data corresponding to pdb file.
- [crystal_name]\_event.(map/ccp4) - This has a higher signal-to-noise ratio which will amplify the evidence of ligand occupancy.
- [crystal_name]\_2fofc.(map/ccp4) - estimate of the true electron density from diffraction data and atomic model.
- [crystal_name]\_fofc.(map/ccp4) - difference electron density map, negative density typically represents where no electron density is found but exists in the atom model. Positive densities represent electron density without mapped atom model.

## extra_files

If this is present the files in this folder will have been added by the uploader of the data and has no defined structure. As a result we cannot guess what the contents of the file may be but we hope that the uploader of the extra files will have provided a similar
Files in this folder will be added by the uploader and are largely freeform. Hopefully there will be a readme inside to describe each of the added files.

## Snapshot Details 
