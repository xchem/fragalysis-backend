# Commands for managing specific tasks in Fragalysis

The commands in this directory can be run from the container directly using
a command of the form:
```
python manage.py <python script> <options>
```

It should be possible to run them either manually or as part of a cron job to 
do maintenance and data intergrity-style tasks.

It would be good if all commands are lightly documented here.

## tags_from_sites

Data correction to create MoleculeTag records from the existing
MolGroup group objects. This is automatically done as part of the 
target_set_upload process, so future records are dealt with. However existing
datasets need to be updated. Doing a full conversion run may have hit 
unexpected errors with old data, so it was agreed to update the tags on a
Target by target basis. 

**Format of command**
```
python manage.py tags_from_sites <target name> <update = yes/no>
```
where:
- update = 'no' will do a dry run outputting the records that will be 
updated if 'yes' was set. Checks the number of records against sites.csv 
- update = 'yes' indicates that the update to the datebase should be performed

**Notes**

1. In case of problems, during testing I removed the tags using the following 
sql commands:
```
begin work;
delete from viewer_moleculetag_molecules vmm 
   where vmm.moleculetag_id in (select id from viewer_moleculetag vm where target_id = <target id>);
delete from viewer_moleculetag vm where target_id = <target id>;
commit;
```

2. For older targets, before the sites functionality was implemented, the 
centre of mass was used as an equivalent grouping. If there is no sites.csv file, the 
function will look for molgroups with a description of 'c_of_m' and create
site tags for those instead with a name 'c_of_m_<molgroup index>'.
