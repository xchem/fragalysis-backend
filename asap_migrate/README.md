# ASAP database migrations

## Create working test database

1. Copy backup from NFS server in prod cluster. Will be found as something like this: 
   /nfs/kubernetes-prod/production-stack-database-backup-pvc-d2affd17-5291-4c19-b4e3-a6253f35a737/hourly/backup-2022-11-28T09:51:11Z-dumpall.sql.gz
2. Gunzip the backup, rename to backup-dumpall.sql and place in this directory.
3. Start the container: `docker-compose up`
4. Shell into the container: `docker exec -it database bash`
5. Restore the backup: `psql -U admin -f /backup-dumpall.sql`
6. Checkout the contents: `psql -U admin -d frag`
7. Checkout the contents as the `fragalysis` user: `psql -U fragalysis -d frag`

## viewer_compound table

```
frag=> \d viewer_compound
                                          Table "public.viewer_compound"
       Column       |          Type          | Collation | Nullable |                   Default                   
--------------------+------------------------+-----------+----------+---------------------------------------------
 id                 | integer                |           | not null | nextval('viewer_compound_id_seq'::regclass)
 inchi              | character varying(255) |           | not null | 
 smiles             | character varying(255) |           | not null | 
 mol_log_p          | double precision       |           | not null | 
 mol_wt             | double precision       |           | not null | 
 tpsa               | double precision       |           | not null | 
 heavy_atom_count   | integer                |           | not null | 
 heavy_atom_mol_wt  | double precision       |           | not null | 
 nhoh_count         | integer                |           | not null | 
 no_count           | integer                |           | not null | 
 num_h_acceptors    | integer                |           | not null | 
 num_h_donors       | integer                |           | not null | 
 num_het_atoms      | integer                |           | not null | 
 num_rot_bonds      | integer                |           | not null | 
 num_val_electrons  | integer                |           | not null | 
 ring_count         | integer                |           | not null | 
 all_identifiers    | text                   |           |          | 
 comments           | text                   |           |          | 
 current_identifier | character varying(255) |           |          | 
 description        | text                   |           |          | 
 long_inchi         | text                   |           |          | 
Indexes:
    "viewer_compound_pkey" PRIMARY KEY, btree (id)
    "viewer_compound_inchi_long_inchi_c05664b2_uniq" UNIQUE CONSTRAINT, btree (inchi, long_inchi)
    "viewer_compound_current_identifier_ce83eec4" btree (current_identifier)
    "viewer_compound_current_identifier_ce83eec4_like" btree (current_identifier varchar_pattern_ops)
    "viewer_compound_inchi_e90306e6" btree (inchi)
    "viewer_compound_inchi_e90306e6_like" btree (inchi varchar_pattern_ops)
    "viewer_compound_smiles_b9cfac6a" btree (smiles)
    "viewer_compound_smiles_b9cfac6a_like" btree (smiles varchar_pattern_ops)
Referenced by:
    TABLE "hypothesis_vector" CONSTRAINT "hypothesis_vector_cmpd_id_id_37cba707_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "scoring_cmpdchoice" CONSTRAINT "scoring_cmpdchoice_cmpd_id_id_e5eebb9e_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_activitypoint" CONSTRAINT "viewer_activitypoint_cmpd_id_id_203ed5b1_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_compound_inspirations" CONSTRAINT "viewer_compound_insp_compound_id_2ae2d1d4_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_compound_project_id" CONSTRAINT "viewer_compound_proj_compound_id_1d7463b1_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_computedmolecule" CONSTRAINT "viewer_computedmolec_compound_id_1b9ae943_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_designset_compounds" CONSTRAINT "viewer_designset_com_compound_id_d2beec89_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_molecule" CONSTRAINT "viewer_molecule_cmpd_id_id_778fdcd7_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
```

Smiles column is duplicated:

```
frag=> select count(id), smiles from viewer_compound group by smiles limit 10;
 count |              smiles              
-------+----------------------------------
     1 | 
     1 | Brc1cccnc1
     1 | Brc1ccc(NN=Nc2ccccc2)cc1
     1 | Brc1ccc(OCCCN2CCOCC2)cc1
    10 | Brc1ccnc2ncccc12
     1 | Brc1cn[nH]c1
     1 | Br[C@@H]1C=CN=C(NCN=CNc2cccs2)C1
     2 | C
     1 | c1cc2c(c(OCc3ccon3)c1)CCCC2
     1 | c1ccc2c(c1)CC1(CCCN1)C2
(10 rows)
```

```
frag=> select count(c.id), c.smiles, min(c.id), array_agg(c.id) from viewer_compound c group by c.smiles limit 20;
 count |              smiles              |  min  |                          array_agg                           
-------+----------------------------------+-------+--------------------------------------------------------------
     1 |                                  |  1609 | {1609}
     1 | Brc1cccnc1                       |  4431 | {4431}
     1 | Brc1ccc(NN=Nc2ccccc2)cc1         | 14762 | {14762}
     1 | Brc1ccc(OCCCN2CCOCC2)cc1         | 17885 | {17885}
    10 | Brc1ccnc2ncccc12                 |  3973 | {3973,10552,10553,10559,10560,10561,10912,15794,16306,16317}
     1 | Brc1cn[nH]c1                     |  1558 | {1558}
     1 | Br[C@@H]1C=CN=C(NCN=CNc2cccs2)C1 |  5830 | {5830}
     2 | C                                |   973 | {973,10077}
     1 | c1cc2c(c(OCc3ccon3)c1)CCCC2      | 18055 | {18055}
     1 | c1ccc2c(c1)CC1(CCCN1)C2          | 16679 | {16679}
     3 | c1cc(-c2ccc3c(c2)OCO3)n[nH]1     |    96 | {96,15762,16656}
     1 | C1CC(C2CCC3OCOC3C2)NN1           |   881 | {881}
     1 | c1ccc2c(CN3CCCCCC3)c[nH]c2c1     |   448 | {448}
     1 | c1cc(C2CCNCC2)cc(C2CNC2)c1       | 14043 | {14043}
     4 | c1cc(C2CCNCC2)n[nH]1             |  6514 | {6514,11048,16565,16566}

```

Create tmp table with the id lookups:
```
frag=> create table tmp_cmpds (id integer, mid integer, smiles varchar(255), c integer);
CREATE TABLE
```

Insert the data:
```
frag=> with t as (select min(id) mid, smiles, count(id) cnt from viewer_compound group by smiles) insert into tmp_cmpds(id, mid, smiles, c) select c.id, t.mid, c.smiles, t.cnt from viewer_compound c join t on t.smiles = c.smiles limit 10;
INSERT 0 10
```

Create index:
```
create index on tmp_cmpds(id);
```

## hypothesis_vector table

```
frag=> \d hypothesis_vector
                                      Table "public.hypothesis_vector"
   Column   |          Type          | Collation | Nullable |                    Default                    
------------+------------------------+-----------+----------+-----------------------------------------------
 id         | integer                |           | not null | nextval('hypothesis_vector_id_seq'::regclass)
 smiles     | character varying(255) |           | not null | 
 type       | character varying(2)   |           | not null | 
 cmpd_id_id | integer                |           | not null | 
Indexes:
    "hypothesis_vector_pkey" PRIMARY KEY, btree (id)
    "hypothesis_vector_cmpd_id_id_smiles_type_e8e2b79c_uniq" UNIQUE CONSTRAINT, btree (cmpd_id_id, smiles, type)
    "hypothesis_vector_cmpd_id_id_37cba707" btree (cmpd_id_id)
Foreign-key constraints:
    "hypothesis_vector_cmpd_id_id_37cba707_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
Referenced by:
    TABLE "hypothesis_vector3d" CONSTRAINT "hypothesis_vector3d_vector_id_id_91d4fe51_fk_hypothesi" FOREIGN KEY (vector_id_id) REFERENCES hypothesis_vector(id) DEFERRABLE INITIALLY DEFERRED
```

```
frag=> select c.id, c.smiles, v.id, v.smiles from viewer_compound c join hypothesis_vector v on v.cmpd_id_id = c.id where c.smiles = 'Brc1ccnc2ncccc12';
  id   |      smiles      |  id   |         smiles         
-------+------------------+-------+------------------------
  3973 | Brc1ccnc2ncccc12 | 48736 | [Xe]c1ccnc2ncccc12
  3973 | Brc1ccnc2ncccc12 | 48737 | Brc1cc([Xe])nc2ncccc12
  3973 | Brc1ccnc2ncccc12 | 48738 | Brc1ccnc2ncc([Xe])cc12
  3973 | Brc1ccnc2ncccc12 | 48739 | Brc1c([Xe])cnc2ncccc12
  3973 | Brc1ccnc2ncccc12 | 48740 | Brc1ccnc2nccc([Xe])c12
  3973 | Brc1ccnc2ncccc12 | 48741 | Brc1ccnc2nc([Xe])ccc12
 10552 | Brc1ccnc2ncccc12 | 61349 | [Xe]c1ccnc2ncccc12
 10552 | Brc1ccnc2ncccc12 | 61350 | Brc1c([Xe])cnc2ncccc12
 10552 | Brc1ccnc2ncccc12 | 61351 | Brc1cc([Xe])nc2ncccc12
 10552 | Brc1ccnc2ncccc12 | 61352 | Brc1ccnc2nc([Xe])ccc12
 10552 | Brc1ccnc2ncccc12 | 61353 | Brc1ccnc2ncc([Xe])cc12
 ...
```

Verification:
```
frag=> select count(*) from hypothesis_vector h join tmp_cmpds t on h.cmpd_id_id = t.id where h.cmpd_id_id != t.mid;
 count 
-------
 17208
(1 row)
```

```
frag=> with t as (select count(id) c from viewer_compound group by smiles) select count(*) from t where c = 1; 
 count 
-------
 16213
(1 row)
```

Update the values:
```
frag=> update hypothesis_vector h set cmpd_id_id = (select t.mid from tmp_cmpds t where h.cmpd_id_id = t.id);
ERROR:  duplicate key value violates unique constraint "hypothesis_vector_cmpd_id_id_smiles_type_e8e2b79c_uniq"
DETAIL:  Key (cmpd_id_id, smiles, type)=(1669, O=C(CS)Nc1nc2c([Xe])cccc2s1, DE) already exists.
```
That table also needs de-duplicating? But how is the data generated? That probably needs updating too?

Proposed action:
1. examine how this data is populated to establish strategy
2. fixing probably involves dropping the constraint, updating IDs, de-duplicating and restoring the constraint


## viewer_compound_project_id

```
frag=> \d viewer_compound_project_id
                               Table "public.viewer_compound_project_id"
   Column    |  Type   | Collation | Nullable |                        Default                         
-------------+---------+-----------+----------+--------------------------------------------------------
 id          | integer |           | not null | nextval('viewer_compound_project_id_id_seq'::regclass)
 compound_id | integer |           | not null | 
 project_id  | integer |           | not null | 
Indexes:
    "viewer_compound_project_id_pkey" PRIMARY KEY, btree (id)
    "viewer_compound_project_id_compound_id_project_id_14bf5e78_uniq" UNIQUE CONSTRAINT, btree (compound_id, project_id)
    "viewer_compound_project_id_compound_id_1d7463b1" btree (compound_id)
    "viewer_compound_project_id_project_id_dd92baa4" btree (project_id)
Foreign-key constraints:
    "viewer_compound_proj_compound_id_1d7463b1_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    "viewer_compound_proj_project_id_dd92baa4_fk_viewer_pr" FOREIGN KEY (project_id) REFERENCES viewer_project(id) DEFERRABLE INITIALLY DEFERRED
```

Proposed action:
1. drop the viewer_compound_project_id_compound_id_project_id_14bf5e78_uniq constraint
2. update the compound_id values to the new ones
3. remove duplicate rows
4. re-appy the unique constraint
5. establish what creates data in this table to ensure correct data is created

## viewer_computedmolecule table

```
frag=> \d viewer_computedmolecule
                                        Table "public.viewer_computedmolecule"
     Column      |          Type          | Collation | Nullable |                       Default                       
-----------------+------------------------+-----------+----------+-----------------------------------------------------
 id              | integer                |           | not null | nextval('viewer_computedmolecule_id_seq'::regclass)
 sdf_info        | text                   |           | not null | 
 name            | character varying(50)  |           | not null | 
 smiles          | character varying(255) |           | not null | 
 compound_id     | integer                |           | not null | 
 computed_set_id | character varying(50)  |           | not null | 
 pdb_id          | integer                |           |          | 
Indexes:
    "viewer_computedmolecule_pkey" PRIMARY KEY, btree (id)
    "viewer_computedmolecule_compound_id_1b9ae943" btree (compound_id)
    "viewer_computedmolecule_computed_set_id_68630526" btree (computed_set_id)
    "viewer_computedmolecule_computed_set_id_68630526_like" btree (computed_set_id varchar_pattern_ops)
    "viewer_computedmolecule_pdb_id_5749b840" btree (pdb_id)
Foreign-key constraints:
    "viewer_computedmolec_compound_id_1b9ae943_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    "viewer_computedmolec_computed_set_id_68630526_fk_viewer_co" FOREIGN KEY (computed_set_id) REFERENCES viewer_computedset(name) DEFERRABLE INITIALLY DEFERRED
    "viewer_computedmolecule_pdb_id_5749b840_fk_viewer_protein_id" FOREIGN KEY (pdb_id) REFERENCES viewer_protein(id) DEFERRABLE INITIALLY DEFERRED
Referenced by:
    TABLE "viewer_computedmolecule_computed_inspirations" CONSTRAINT "viewer_computedmolec_computedmolecule_id_2398f8f5_fk_viewer_co" FOREIGN KEY (computedmolecule_id) REFERENCES viewer_computedmolecule(id) DEFERRABLE INITIALL
Y DEFERRED
    TABLE "viewer_numericalscorevalues" CONSTRAINT "viewer_numericalscor_compound_id_2b2593f3_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_computedmolecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_textscorevalues" CONSTRAINT "viewer_textscorevalu_compound_id_ecb14f5a_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_computedmolecule(id) DEFERRABLE INITIALLY DEFERRED
```

No unique constraints. Process seems straight forward:
```
frag=> update viewer_computedmolecule h set compound_id = (select t.mid from tmp_cmpds t where h.compound_id = t.id);
UPDATE 5087
```

Proposed action:
1. establish how this table is populated (computed set upload?) to ensure that correct compound_id values are used. 


## scoring_cmpdchoice table

```
frag=> \d scoring_cmpdchoice
Table "public.scoring_cmpdchoice"
Column    |         Type         | Collation | Nullable |                    Default                     
-------------+----------------------+-----------+----------+------------------------------------------------
id          | integer              |           | not null | nextval('scoring_cmpdchoice_id_seq'::regclass)
choice_type | character varying(2) |           | not null |
score       | double precision     |           |          |
cmpd_id_id  | integer              |           | not null |
user_id_id  | integer              |           | not null |
Indexes:
"scoring_cmpdchoice_pkey" PRIMARY KEY, btree (id)
"scoring_cmpdchoice_user_id_id_cmpd_id_id_ch_f95fbaaa_uniq" UNIQUE CONSTRAINT, btree (user_id_id, cmpd_id_id, choice_type)
"scoring_cmpdchoice_cmpd_id_id_e5eebb9e" btree (cmpd_id_id)
"scoring_cmpdchoice_user_id_id_82d4600d" btree (user_id_id)
Foreign-key constraints:
"scoring_cmpdchoice_cmpd_id_id_e5eebb9e_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
"scoring_cmpdchoice_user_id_id_82d4600d_fk_auth_user_id" FOREIGN KEY (user_id_id) REFERENCES auth_user(id) DEFERRABLE INITIALLY DEFERRED
```

Currently contains no data so no data updates needed.

Proposed action:
1. establish the purpose of this table. If it is used in the Fraglysis code then some action might be needed.

## viewer_activitypoint table

```
frag=> \d viewer_activitypoint
                                       Table "public.viewer_activitypoint"
    Column    |          Type          | Collation | Nullable |                     Default                      
--------------+------------------------+-----------+----------+--------------------------------------------------
 id           | integer                |           | not null | nextval('viewer_activitypoint_id_seq'::regclass)
 source       | character varying(50)  |           |          | 
 activity     | double precision       |           | not null | 
 units        | character varying(50)  |           | not null | 
 confidence   | integer                |           |          | 
 internal_id  | character varying(150) |           |          | 
 operator     | character varying(5)   |           | not null | 
 cmpd_id_id   | integer                |           | not null | 
 target_id_id | integer                |           | not null | 
Indexes:
    "viewer_activitypoint_pkey" PRIMARY KEY, btree (id)
    "viewer_activitypoint_target_id_id_activity_cm_11d84f70_uniq" UNIQUE CONSTRAINT, btree (target_id_id, activity, cmpd_id_id, units)
    "viewer_activitypoint_activity_c36634de" btree (activity)
    "viewer_activitypoint_cmpd_id_id_203ed5b1" btree (cmpd_id_id)
    "viewer_activitypoint_confidence_4a5bafc8" btree (confidence)
    "viewer_activitypoint_source_2e56854c" btree (source)
    "viewer_activitypoint_source_2e56854c_like" btree (source varchar_pattern_ops)
```

Currently contains no data but if it ever is used the unique constraint needs to be investigated.

Proposed action:
1. establish the purpose of this table. If it is used in the Fraglysis code then some action might be needed.

## viewer_compound_inspirations table

```
frag=> \d viewer_compound_inspirations
                               Table "public.viewer_compound_inspirations"
   Column    |  Type   | Collation | Nullable |                         Default                          
-------------+---------+-----------+----------+----------------------------------------------------------
 id          | integer |           | not null | nextval('viewer_compound_inspirations_id_seq'::regclass)
 compound_id | integer |           | not null | 
 molecule_id | integer |           | not null | 
Indexes:
    "viewer_compound_inspirations_pkey" PRIMARY KEY, btree (id)
    "viewer_compound_inspirat_compound_id_molecule_id_b02d4298_uniq" UNIQUE CONSTRAINT, btree (compound_id, molecule_id)
    "viewer_compound_inspirations_compound_id_2ae2d1d4" btree (compound_id)
    "viewer_compound_inspirations_molecule_id_023fdb35" btree (molecule_id)
Foreign-key constraints:
    "viewer_compound_insp_compound_id_2ae2d1d4_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    "viewer_compound_insp_molecule_id_023fdb35_fk_viewer_mo" FOREIGN KEY (molecule_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
```

Currently contains no data.

Proposed action:
1. establish the purpose of this table. If it is used in the Fraglysis code then some action might be needed.

## viewer_designset_compounds table

```
frag=> \d viewer_designset_compounds
                               Table "public.viewer_designset_compounds"
    Column    |  Type   | Collation | Nullable |                        Default                         
--------------+---------+-----------+----------+--------------------------------------------------------
 id           | integer |           | not null | nextval('viewer_designset_compounds_id_seq'::regclass)
 designset_id | integer |           | not null | 
 compound_id  | integer |           | not null | 
Indexes:
    "viewer_designset_compounds_pkey" PRIMARY KEY, btree (id)
    "viewer_designset_compoun_designset_id_compound_id_d84f733d_uniq" UNIQUE CONSTRAINT, btree (designset_id, compound_id)
    "viewer_designset_compounds_compound_id_d2beec89" btree (compound_id)
    "viewer_designset_compounds_designset_id_0861bc82" btree (designset_id)
Foreign-key constraints:
    "viewer_designset_com_compound_id_d2beec89_fk_viewer_co" FOREIGN KEY (compound_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    "viewer_designset_com_designset_id_0861bc82_fk_viewer_de" FOREIGN KEY (designset_id) REFERENCES viewer_designset(id) DEFERRABLE INITIALLY DEFERRED
```

Currently contains no data.

Proposed action:
1. establish the purpose of this table. If it is used in the Fraglysis code then some action might be needed.

## viewer_molecule table

```
frag=> \d viewer_molecule
                                      Table "public.viewer_molecule"
   Column   |          Type          | Collation | Nullable |                   Default                   
------------+------------------------+-----------+----------+---------------------------------------------
 id         | integer                |           | not null | nextval('viewer_molecule_id_seq'::regclass)
 smiles     | character varying(255) |           |          | 
 lig_id     | character varying(5)   |           |          | 
 chain_id   | character varying(1)   |           |          | 
 mol_type   | character varying(2)   |           | not null | 
 sdf_info   | text                   |           |          | 
 rscc       | double precision       |           |          | 
 occupancy  | double precision       |           |          | 
 x_com      | double precision       |           |          | 
 y_com      | double precision       |           |          | 
 z_com      | double precision       |           |          | 
 rmsd       | double precision       |           |          | 
 cmpd_id_id | integer                |           | not null | 
 prot_id_id | integer                |           | not null | 
 sdf_file   | character varying(255) |           |          | 
Indexes:
    "viewer_molecule_pkey" PRIMARY KEY, btree (id)
    "viewer_molecule_prot_id_id_cmpd_id_id_mol_type_6e3feda6_uniq" UNIQUE CONSTRAINT, btree (prot_id_id, cmpd_id_id, mol_type)
    "viewer_molecule_cmpd_id_id_778fdcd7" btree (cmpd_id_id)
    "viewer_molecule_prot_id_id_26e598bd" btree (prot_id_id)
    "viewer_molecule_smiles_20a42f4c" btree (smiles)
    "viewer_molecule_smiles_20a42f4c_like" btree (smiles varchar_pattern_ops)
Foreign-key constraints:
    "viewer_molecule_cmpd_id_id_778fdcd7_fk_viewer_compound_id" FOREIGN KEY (cmpd_id_id) REFERENCES viewer_compound(id) DEFERRABLE INITIALLY DEFERRED
    "viewer_molecule_prot_id_id_26e598bd_fk_viewer_protein_id" FOREIGN KEY (prot_id_id) REFERENCES viewer_protein(id) DEFERRABLE INITIALLY DEFERRED
Referenced by:
    TABLE "hypothesis_interactionpoint" CONSTRAINT "hypothesis_interacti_mol_id_id_4d4038fb_fk_viewer_mo" FOREIGN KEY (mol_id_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "hypothesis_vector3d" CONSTRAINT "hypothesis_vector3d_mol_id_id_810911a3_fk_viewer_molecule_id" FOREIGN KEY (mol_id_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "scoring_molannotation" CONSTRAINT "scoring_molannotation_mol_id_id_7a72d48f_fk_viewer_molecule_id" FOREIGN KEY (mol_id_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "scoring_molchoice" CONSTRAINT "scoring_molchoice_mol_id_id_cb0599f0_fk_viewer_molecule_id" FOREIGN KEY (mol_id_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "scoring_molgroup_mol_id" CONSTRAINT "scoring_molgroup_mol_molecule_id_546cf5c3_fk_viewer_mo" FOREIGN KEY (molecule_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "scoring_scorechoice" CONSTRAINT "scoring_scorechoice_mol_id_id_dd2ad6a7_fk_viewer_molecule_id" FOREIGN KEY (mol_id_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_compound_inspirations" CONSTRAINT "viewer_compound_insp_molecule_id_023fdb35_fk_viewer_mo" FOREIGN KEY (molecule_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_computedmolecule_computed_inspirations" CONSTRAINT "viewer_computedmolec_molecule_id_57c97346_fk_viewer_mo" FOREIGN KEY (molecule_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
    TABLE "viewer_moleculetag_molecules" CONSTRAINT "viewer_moleculetag_m_molecule_id_42340bc9_fk_viewer_mo" FOREIGN KEY (molecule_id) REFERENCES viewer_molecule(id) DEFERRABLE INITIALLY DEFERRED
```

Update:
```
frag=> update viewer_molecule m set cmpd_id_id = (select t.mid from tmp_cmpds t where m.cmpd_id_id = t.id);
UPDATE 4320
```

Proposed action:
1. investigate how this table is populated (target loader?) to ensure that compound refs are unique.