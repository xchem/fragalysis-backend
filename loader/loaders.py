import os,sys,re,json
import shutil
from viewer.models import Target,Protein,Molecule,Compound
from pandda.models import PanddaSite,PanddaEvent
from hypothesis.models import Vector3D,Vector,InteractionPoint,TargetResidue,ProteinResidue,Interaction
from hypothesis.definitions import VectTypes,IntTypes
from django.core.exceptions import ValidationError
from django.core.files import File
from rdkit import Chem
from rdkit.Chem import Descriptors
from scoring.models import MolGroup
from frag.alysis.run_clustering import run_lig_cluster
from loader.functions import sanitize_mol,get_path_or_none
from frag.network.decorate import get_3d_vects_for_mol

def add_target(title):
    """
    Add a target
    :param title: add a target by titlte
    :return: the created target
    """
    new_target = Target.objects.get_or_create(title=title)[0]
    return new_target

def add_prot(pdb_file_path, code, target, mtz_path=None, map_path=None):
    """
    Add a protein with a PDB, code and
    :param pdb_file_path: the PDB file path
    :param code: the unique code for this file
    :param target: the target to be linkede to
    :param mtz_path: the path to the MTZ file
    :param map_path: the path to the MAP file
    :return: the created protein
    """
    new_prot = Protein.objects.get_or_create(code=code,target_id=target)[0]
    new_prot.apo_holo = True
    new_prot.pdb_info.save(os.path.basename(pdb_file_path), File(open(pdb_file_path)))
    if mtz_path:
        new_prot.mtz_info.save(os.path.basename(mtz_path),File(open(mtz_path)))
    if map_path:
        new_prot.map_info.save(os.path.basename(map_path),File(open(map_path)))
    new_prot.save()
    return new_prot

def add_comp(mol, option=None, comp_id=None):
    """
    Function to add a new compound to the database given an RDKit molecule
    Takes an RDKit molecule. Option of LIG to return original smiles with the Compound object
    Returns a compound object for the RDKit molecule.
    :param mol:
    :param option:
    :param comp_id:
    :return:
    """
    # Neutralise and desalt compound the compound
    sanitized_mol = sanitize_mol(mol)
    # Store the isomeric smiles
    smiles = Chem.MolToSmiles(sanitized_mol, isomericSmiles=True)
    # The inchi string is used for unique identification
    inchi = Chem.MolToInchi(sanitized_mol)
    # Now convert back to inchi to canonicalise
    tmp_mol = Chem.MolFromInchi(inchi)
    if tmp_mol is None:
        # If error in INNCHI READ -> NOT NECCESARILY A KILLER
        sys.stderr.write("INCHI ERROR: " + inchi)
    else:
        inchi = Chem.MolToInchi(tmp_mol)
    # Now attribute all this meta-deta to the compound object
    new_comp = Compound()
    new_comp.smiles = smiles
    if len(smiles) > Compound._meta.get_field('smiles').max_length:
        print "SMILES TOO LONG"
        return None
    new_comp.inchi = inchi
    m = sanitized_mol
    try:
        new_comp.mol_log_p = Chem.Crippen.MolLogP(m)
    except:
        sys.stderr.write("NONE MOLECULE PRODUCED\n" + smiles + "\n" + inchi)
        return None
    new_comp.mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(m))
    new_comp.heavy_atom_count = Chem.Lipinski.HeavyAtomCount(m)
    new_comp.heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(m))
    new_comp.nhoh_count = Chem.Lipinski.NHOHCount(m)
    new_comp.no_count = Chem.Lipinski.NOCount(m)
    new_comp.num_h_acceptors = Chem.Lipinski.NumHAcceptors(m)
    new_comp.num_h_donors = Chem.Lipinski.NumHDonors(m)
    new_comp.num_het_atoms = Chem.Lipinski.NumHeteroatoms(m)
    new_comp.num_rot_bonds = Chem.Lipinski.NumRotatableBonds(m)
    new_comp.num_val_electrons = Descriptors.NumValenceElectrons(m)
    new_comp.ring_count = Chem.Lipinski.RingCount(m)
    new_comp.tpsa = Chem.rdMolDescriptors.CalcTPSA(m)
    # Validate that the compound is unique
    try:
        new_comp.validate_unique()
        new_comp.save()
        return new_comp
    except ValidationError:
        return Compound.objects.get(inchi=inchi)

def add_mol(mol_sd,prot,lig_id="LIG",chaind_id="Z",occupancy=0.0):
    """
    Function to add a new Molecule to the database
    :param mol_sd: the SDMolBlock of the molecule
    :param prot: the protein it is associated to
    :param lig_id: the 3 letter ligand id
    :param chaind_id: the chain id
    :param occupancy: the occupancy
    :return: the created molecule
    """
    rd_mol = Chem.MolFromMolFile(mol_sd)
    if rd_mol is None:
        return None
    # Get the reference compound
    comp_ref = add_comp(rd_mol)
    new_mol = Molecule.objects.get_or_create(prot_id=prot,cmpd_id=comp_ref)[0]
    # Make a protein object by which it is related in the DB
    new_mol.sdf_info = Chem.MolToMolBlock(rd_mol)
    new_mol.smiles = Chem.MolToSmiles(rd_mol, isomericSmiles=True)
    # Find out how to add this information from Proasis
    new_mol.lig_id = lig_id
    new_mol.chain_id = chaind_id
    new_mol.occupancy = occupancy
    # Add this to the compound list -> make sure this passes in for the
    # correct molecule. I.e. if it fails where does it go???
    # Now link that compound back
    new_mol.cmpd_id = comp_ref
    new_mol.save()
    return new_mol


def parse_proasis(input_string):
    return input_string[:3].strip(), int(input_string[5:].strip()), input_string[3:5].strip()

def add_contacts(input_data,target,prot,mol):
    int_type = IntTypes()
    if type(input_data)== dict:
        int_list = input_data["results"]
    else:
        int_list = input_data
    for interaction in int_list:
        # Ignore Water mediasted hypothesis with Protein for now
        if interaction['hetmoltype'] == 'WATER':
            continue
        res_name, res_num, chain_id = parse_proasis(interaction['dstrname'])
        targ_res = TargetResidue.objects.get_or_create(target_id=target,res_name=res_name,res_num=res_num,chain_id=chain_id)[0]
        prot_res = ProteinResidue.objects.get_or_create(targ_res_id=targ_res,prot_id=prot)[0]
        interation_point = InteractionPoint.objects.get_or_create(prot_res_id=prot_res,mol_id=mol,
                                                                  protein_atom_name=interaction['dstaname'],
                                                                  molecule_atom_name=interaction['srcaname']
                                                                  )[0]
        Interaction.objects.get_or_create(interaction_version="PR",interaction_point=interation_point,
                                          interaction_type=int_type.get_int_conv("PR",interaction["contactType"]),
                                          distance=interaction["dis"],prot_smarts=interaction['dstType'],
                                          mol_smarts=interaction['srcType'])

def load_from_dir(target_name, dir_path, input_dict):
    """
    Load the data for a given target from a directory structure
    :param target_name:
    :param dir_path:
    :return:
    """
    if(os.path.isdir(dir_path)):
        pass
    else:
        "No data to add: " + target_name
        return None
    new_target = add_target(target_name)
    directories = os.listdir(dir_path)
    for xtal in directories:
        print(xtal)
        new_path = os.path.join(dir_path, xtal)
        pdb_file_path = get_path_or_none(new_path,xtal,input_dict,"APO")
        mol_file_path = get_path_or_none(new_path,xtal,input_dict,"MOL")
        map_path = get_path_or_none(new_path,xtal,input_dict,"MAP")
        mtz_path = get_path_or_none(new_path,xtal,input_dict,"MTZ")
        contact_path = get_path_or_none(new_path,xtal,input_dict,"CONTACTS")
        if not pdb_file_path or not mol_file_path:
            continue
        if os.path.isfile(pdb_file_path) and os.path.isfile(mol_file_path):
            new_prot = add_prot(pdb_file_path,xtal,new_target,mtz_path=mtz_path,map_path=map_path)
            new_mol = add_mol(mol_file_path, new_prot)
            if not new_mol:
                print("NONE MOL: "+xtal)
            else:
                if contact_path:
                    if os.path.isfile(contact_path):
                        add_contacts(json.load(open(contact_path)),new_target,new_prot,new_mol)
        else:
            print("File not found: "+xtal)


def parse_centre(input_str):
    return json.loads(input_str.strip('"'))


def create_site(site,target,pandda_version,site_align_cent,site_native_cent):
    new_site = PanddaSite.objects.get_or_create(site_id=site, target_id=target,pandda_run="STANDARD")[0]
    new_site.pandda_version = pandda_version
    new_site.site_align_com_x = site_align_cent[0]
    new_site.site_align_com_y = site_align_cent[1]
    new_site.site_align_com_z = site_align_cent[2]
    new_site.site_native_com_x = site_native_cent[0]
    new_site.site_native_com_y = site_native_cent[1]
    new_site.site_native_com_z = site_native_cent[2]
    new_site.save()
    return new_site

def create_event(xtal,event,site,pandda_version,pdb_file,mtz_path,map_path,lig_id,
                     event_cent,event_dist,lig_cent,lig_dist,site_align_cent,site_native_cent,target):
    # Now make the event
    new_site = create_site(site,target,pandda_version,site_align_cent,site_native_cent)
    new_event = PanddaEvent.objects.get_or_create(xtal=xtal, event=event, pandda_site=new_site, target_id=target)[0]
    new_event.pdb_info.save(os.path.basename(pdb_file), File(open(pdb_file)))
    new_event.mtz_info.save(os.path.basename(mtz_path),File(open(mtz_path)))
    new_event.map_info.save(os.path.basename(map_path),File(open(map_path)))
    small_map_path = map_path.replace(".map",'_small.map')
    new_event.small_map_info.save(os.path.basename(small_map_path),File(open(small_map_path)))
    new_event.lig_id = lig_id
    new_event.event_com_x = event_cent[0]
    new_event.event_com_y = event_cent[1]
    new_event.event_com_z = event_cent[2]
    new_event.lig_com_x = lig_cent[0]
    new_event.lig_com_y = lig_cent[1]
    new_event.lig_com_z = lig_cent[2]
    new_event.event_dist_from_site_centroid = event_dist
    new_event.lig_dist_from_site_centroid = lig_dist
    new_event.save()
    return new_event,new_site

def load_events_from_dir(target_name,dir_path):
    new_target = add_target(target_name)
    if(os.path.isdir(dir_path)):
        os.chdir(dir_path)
    else:
        "No events to add: " + target_name
        return None
    lines = [x.strip() for x in open("out.csv").readlines()]
    header = lines[0].split(",")
    PATTERN = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    for line in lines[1:]:
        spl_line = PATTERN.split(line)[1::2]
        xtal = spl_line[0]
        event = spl_line[1]
        site = spl_line[2]
        pandda_version = spl_line[3]
        pdb_file = spl_line[4]
        mtz_file = spl_line[5]
        map_file = spl_line[6]
        lig_id = spl_line[7]
        event_cent = parse_centre(spl_line[8])
        event_dist = float(spl_line[9])
        lig_cent = parse_centre(spl_line[10])
        lig_dist = float(spl_line[11])
        site_align_cent = parse_centre(spl_line[12])
        site_native_cent = parse_centre(spl_line[13])
        create_event(xtal,event,site,pandda_version,pdb_file,mtz_file,map_file,lig_id,
                     event_cent,event_dist,lig_cent,lig_dist,site_align_cent,site_native_cent,new_target)
    os.chdir("../")

def get_vectors(mols):
    """
    Get the vectors for a given molecule
    :param mols:
    :param target:
    :return:
    """
    vect_types = VectTypes()
    for mol in mols:
        if "." in mol.smiles:
            print("SKIPPING - FRAGMENT: " + str(mol.pk)) + " " + str(mol.smiles)
            continue
        vectors = get_3d_vects_for_mol(mol.sdf_info)
        for vect_type in vectors:
            vect_choice = vect_types.translate_vect_types(vect_type)
            for vector in vectors[vect_type]:
                spl_vect = vector.split("__")
                smiles = spl_vect[0]
                if len(spl_vect) > 1:
                    vect_ind = int(spl_vect[1])
                else:
                    vect_ind = 0
                new_vect = Vector.objects.get_or_create(smiles=smiles,cmpd_id=mol.cmpd_id,type=vect_choice)[0]
                new_vect3d = Vector3D.objects.get_or_create(mol_id=mol,vector_id=new_vect,number=vect_ind)[0]
                # The start position
                new_vect3d.start_x = float(vectors[vect_type][vector][0][0])
                new_vect3d.start_y = float(vectors[vect_type][vector][0][1])
                new_vect3d.start_z = float(vectors[vect_type][vector][0][2])
                # The end position
                new_vect3d.end_x = float(vectors[vect_type][vector][1][0])
                new_vect3d.end_y = float(vectors[vect_type][vector][1][1])
                new_vect3d.end_z = float(vectors[vect_type][vector][1][2])
                new_vect3d.save()

def cluster_mols(rd_mols,mols,target):
    id_mols = [x.pk for x in mols]
    out_data = run_lig_cluster(rd_mols, id_mols)
    for clust_type in out_data:
        for cluster in out_data[clust_type]:
            mol_group = MolGroup()
            if clust_type!="c_of_m":
                mol_group.group_type = "PC"
            else:
                mol_group.group_type = "MC"
            mol_group.target_id = target
            mol_group.x_com = out_data[clust_type][cluster]['centre_of_mass'][0]
            mol_group.y_com = out_data[clust_type][cluster]['centre_of_mass'][1]
            mol_group.z_com = out_data[clust_type][cluster]['centre_of_mass'][2]
            mol_group.description = clust_type
            mol_group.save()
            for mol_id in out_data[clust_type][cluster]["mol_ids"]:
                this_mol = Molecule.objects.get(id=mol_id)
                mol_group.mol_id.add(this_mol)

def analyse_mols(mols, target):
    """
    Analyse a list of molecules for a given target
    :param mols:
    :param target:
    :return:
    """
    rd_mols = [Chem.MolFromMolBlock(x.sdf_info) for x in mols]
    cluster_mols(rd_mols,mols,target)
    get_vectors(mols)

def analyse_target(target_name):
    """
    Analyse all the molecules for a particular target
    :param target_name: the name of the target
    :return:
    """
    target = Target.objects.get(title=target_name)
    mols = list(Molecule.objects.filter(prot_id__target_id=target))
    # Delete the old ones for this target
    MolGroup.objects.filter(group_type="PC",target_id=target).delete()
    MolGroup.objects.filter(group_type="MC",target_id=target).delete()
    analyse_mols(mols=mols,target=target)


import csv,os,shutil

FILE_PATH_DICT = {"APO": "_apo.pdb", "MOL": ".mol", "EVENT": "_event.map",
                  "MTZ": ".mtz", "CONTACTS": '_contacts.json'}

def prepare_from_csv(file_path):
    date = "20180430"
    if not os.path.isdir(date):
        os.mkdir(date)
    rows = csv.DictReader(open(file_path))
    for row in rows:
        target_name = row['apo_name'].split("-")[0]
        target_base = os.path.join(date,target_name)
        xtal_name = row['apo_name'].split("_apo.pdb")[0]
        # Make the target directory
        if not os.path.isdir(target_base):
            os.mkdir(target_base)
        xtal_base = os.path.join(target_base,xtal_name)
        # Make this dataset directory
        if not os.path.isdir(xtal_base):
            os.mkdir(xtal_base)
        # Now copy the files over
        shutil.copyfile(os.path.join(row['root_dir'],row['apo_name']),
                        os.path.join(xtal_base,xtal_name+FILE_PATH_DICT["APO"]))
        shutil.copyfile(os.path.join(row['root_dir'], row['mol_name']),
                        os.path.join(xtal_base, xtal_name + FILE_PATH_DICT["MOL"]))
        shutil.copyfile(os.path.join(row['root_dir'], row['fofc_name']),
                        os.path.join(xtal_base, xtal_name + FILE_PATH_DICT["EVENT"]))
        shutil.copyfile(os.path.join(row['root_dir'], row['apo_name'].replace("_apo.pdb","_interactions.json")),
                        os.path.join(xtal_base, xtal_name + FILE_PATH_DICT["CONTACTS"]))
        if row['mtz_name']:
            shutil.copyfile(os.path.join(row['root_dir'], row['mtz_name']),
                            os.path.join(xtal_base, xtal_name + FILE_PATH_DICT["MTZ"]))
# Use this and add into LUIGI pipeline
# prepare_from_csv("/dls/science/groups/i04-1/software/luigi_pipeline/pipeline/logs/proasis_out/proasis_out_20180430.csv")


def process_target(prefix,target_name):
    load_from_dir(target_name, prefix + target_name, FILE_PATH_DICT)
    load_events_from_dir(target_name, prefix + target_name + "/events")
    analyse_target(target_name)