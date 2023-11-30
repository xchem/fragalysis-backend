# Download Readme Template

The download_readme_template file is copied and used as a template for the Directory structure section of the documentation included in the zip files created by the download structures API:
*api/download_structures.

It is placed in an overall structure as follows:

# Documentation for the downloaded zipfile  [Top heading]
## Download details  [subheading]
### Download URLs  [sub-subheading]
(** Note to be updated**)
_Download URL_: https://fragalysis.xchem.diamond.ac.uk/viewer/react/download/tag/a7ea4b13-90b2-4040-8396-6d3fe7b111a3
(** Note to be added**)
_Download snapshot_:  https://https://fragalysis.xchem.diamond.ac.uk/viewer/react/projects/1350/1010"
`[Ensure they render as clickable hyperlink in the pdf, and include https://]`

### Download options selected
(** Note to be added**)
The following options were checked in the download dialogue:
* Selected: All structures
* Selected: PanDDA Event maps - primary evidence
* Selected: Incremental - always up-to-date with latest structures
* Selected: Single SDF of all ligands

[use actual words in the modal - just type them across if necessary]

### Download command (JSON)  
JSON command sent from front-end to backend to generate the download.  This can be reused programmatically as a POST command

The actual JSON (**Note: To be updated currently taken from request, but should be taken from front end**)

## Directory structure  [subheading]
`[The text from download_readme_template.md.]`

## Files included
`[Bulleted List of files in zip]`



- The Download URL
- The API parameters used
- A list of the files in the zip file.

The file will then be turned into a PDF before being added to the zip file.
