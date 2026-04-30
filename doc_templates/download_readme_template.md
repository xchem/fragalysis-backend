## Directory structure

A Fragalysis download will contain a minimum of 2 directories, `aligned_files` and `crytallographic_files`. The download will typically include the additional directories `extra_files`, `scripts` and `yaml_files`, as well as some additional files at the top level directory.

Two important top level files are `metadata.csv` and `smiles.smi`. These are both plain-text files. `metadata.csv` will contain information about the context of each ligand and may provide a convenient way to browse through smiles, site labels and PDB codes for each ligand. `smiles.smi` contains a list of all smiles strings that you have downloaded separated by commans. `[target-name]_combined.sdf` may also be present which will contain all the ligand sdf files in a single sdf file.

### Aligned directory

The aligned directory contains a subdirectory for each dataset that was selected for downloading, aligned to a common reference through [XChemAlign](https://xchem-align.readthedocs.io) processing as they appear in the viewer interface. Depending on your selection of options when downloading the data, the follow file suffixes may be present:

- `[target-name_crystal-name].pdb --- Full atomic model. Protein, ligand, and water/ion/buffer molecules
- `[target-name_crystal-name]_delig-desolv.pdb --- Protein model only. Ligand and water/ion/buffer molecules removed
- `[target-name_crystal-name]_delig-solv.pdb --- Water/ion/buffer molecules molecules only. Protein and ligand molecules removed
- `[target-name_crystal-name]_delig.pdb --- Protein and solvent/ion/buffer molecules. Ligand molecules removed
- `[target-name_crystal-name]_event.ccp4 --- PanDDA event electron density map cut to around 12 Å around the ligand - Background-corrected reflection data higher signal-to-noise enhances ligand evidence corresponding to the PDB file
- `[target-name_crystal-name]_sigmaa.ccp4 --- 2mFo-DFc σA-weighted map cut to around 12 Å around the ligand - Estimate of the true electron density from diffraction data and atomic model
- `[target-name_crystal-name]_diff.ccp4 --- mFo-DFc σA-weighted difference map cut to around 12 Å around the ligand - Negative density indicates model without supporting density, positive density indicates unmodelled features
- `[target-name_crystal-name]_event_crystallographic.ccp4 --- PanDDA event electron density map cut to around 12 Å around the ligand - Background-corrected reflection data higher signal-to-noise enhances ligand evidence corresponding to the PDB file
- `[target-name_crystal-name]_sigmaa_crystallographic.ccp4 --- 2mFo-DFc σA-weighted map cut to around 12 Å around the ligand - Estimate of the true electron density from diffraction data and atomic model
- `[target-name_crystal-name]_diff_crystallographic.ccp4 --- mFo-DFc σA-weighted difference map cut to around 12 Å around the ligand - Negative density indicates model without supporting density, positive density indicates unmodelled features
- `[target-name_crystal-name]_ligand.pdb --- Ligand structure in PDB format
- `[target-name_crystal-name]_ligand.sdf --- Ligand structure in SDF format
- `[target-name_crystal-name]_ligand.smi --- Ligand structure in SMILES format

**IMPORTANT**  
`.ccp4` maps are optimised to work with NGL viewer.  
If viewing in PyMOL or COOT, files that align with the XCA aligned model have the suffix `_crystallographic.ccp4`.

### Crystallographic directory

The `crystallographic_files` directory contains versions of data found in the aligned folder prior to [XChemAlign](https://xchem-align.readthedocs.io) processing. Depending on your selection of options when downloading the data the follow file suffixes may be present:

- `[crystal-name].pdb` --- Full atomic model. Protein, ligand(s), and solvent/ion/buffer molecules                                                                                                         |
- `[crystal-name].mtz` --- Reflection data corresponding to the PDB file                                                                                                |
- `[crystal-name].cif` --- Ligand structure in CIF format                                                                                               |


### Extra files

If the SoakDB CSV and/or SQLite option(s) have been selected, their corresponding files can be found in this directory. Beyond this, if this directory is present the files will have been added by the uploader of the data, and therefore has no defined structure. As a result we cannot guess what the contents of the file may be, but we hope that the uploader of the extra files will have provided a readme inside to describe each of the added files. Some examples of extra files:

- `protein-sequence.fasta`         | Target sequence in FASTA format                                                                                                        |
- `soakdb_[session_number].sqlite` | SoakDB file in SQLite format - Experimental details for each crystal, including soaking conditions, data collection parameters, and processing results. |
- `soakdb_[session_number].csv`    | SoakDB file in CSV format - Experimental details for each crystal, including soaking conditions, data collection parameters, and processing results.           |