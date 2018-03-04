import os,sys
from viewer.models import Target,Protein,Molecule,Compound
from django.core.exceptions import ValidationError
from django.core.files import File
from rdkit import Chem
from rdkit.Chem import Descriptors
from scoring.models import MolGroup
from frag.alysis.run_clustering import run_lig_cluster
from loader.functions import sanitize_mol,get_path_or_none

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

def load_from_dir(target_name, dir_path):
    """
    Load the data for a given target from a directory structure
    :param target_name:
    :param dir_path:
    :return:
    """
    new_target = add_target(target_name)
    directories = os.listdir(dir_path)
    for xtal in directories:
        new_path = os.path.join(dir_path, xtal)
        pdb_file_path = os.path.join(new_path,xtal+"_apo.pdb")
        mol_file_path = os.path.join(new_path,xtal+".mol")
        map_path = get_path_or_none(new_path,xtal,"_event.map")
        mtz_path = get_path_or_none(new_path,xtal,".mtz")
        if os.path.isfile(pdb_file_path) and os.path.isfile(mol_file_path):
            new_prot = add_prot(pdb_file_path,xtal,new_target,mtz_path=mtz_path,map_path=map_path)
            new_mol = add_mol(mol_file_path, new_prot)
            if not new_mol:
                print("NONE MOL: "+xtal)
        else:
            print("File not found: "+xtal)

def analyse_mols(mols, target):
    """
    Analyse a list of molecules for a given target
    :param mols:
    :param target:
    :return:
    """
    rd_mols = [Chem.MolFromMolBlock(x.sdf_info) for x in mols]
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
            mol_group.save()
            for mol_id in out_data[clust_type][cluster]["mol_ids"]:
                this_mol = Molecule.objects.get(id=mol_id)
                mol_group.mol_id.add(this_mol)


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
    analyse_mols(mols=mols)