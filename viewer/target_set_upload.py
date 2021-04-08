"""
target_setupload.py Functions for uploading a Target Dataset to Fragalysis.

This was originally sourced from the code in the fragment-loader repo:
loaders.py
functions.py
"""
import sys, json, os, glob, shutil
from django.contrib.auth.models import User
from viewer.models import Target, Protein, Molecule, Compound, Project, ComputedMolecule
from hypothesis.models import (
    Vector3D,
    Vector,
    InteractionPoint,
    TargetResidue,
    ProteinResidue,
    Interaction,
)
from hypothesis.definitions import VectTypes, IntTypes
from hotspots.models import HotspotMap
from django.core.exceptions import ValidationError
from django.core.files import File
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski, AllChem
from scoring.models import MolGroup,MolAnnotation
from frag.alysis.run_clustering import run_lig_cluster
from frag.network.decorate import get_3d_vects_for_mol
from loader.config import get_dict
import numpy as np
import pandas as pd

from django.conf import settings
from django.core.files.storage import default_storage


# Contribution to the RDKit from Hans de Winter
def _InitialiseNeutralisationReactions():
    """Contribution from Hans de Winter"""
    patts = (
        # Imidazoles
        ("[n+;H]", "n"),
        # Amines
        ("[N+;!H0]", "N"),
        # Carboxylic acids and alcohols
        ("[$([O-]);!$([O-][#7])]", "O"),
        # Thiols
        ("[S-;X1]", "S"),
        # Sulfonamides
        ("[$([N-;X2]S(=O)=O)]", "N"),
        # Enamines
        ("[$([N-;X2][C,N]=C)]", "N"),
        # Tetrazoles
        ("[n-]", "[nH]"),
        # Sulfoxides
        ("[$([S-]=O)]", "S"),
        # Amides
        ("[$([N-]C=O)]", "N"),
    )
    return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in patts]


_reactions = None


def neutralise_charges(smiles, reactions=None):
    """Contribution from Hans de Winter"""
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions = _InitialiseNeutralisationReactions()
        reactions = _reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i, (reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return Chem.MolToSmiles(mol, True), True
    else:
        return smiles, False


def desalt_compound(smiles):
    """Function to desalt compound a given smiles string

    Takes a smiles string
    Returns a desalted smiles string.
    """
    # Chose the biggest fragment, after splitting into fragments
    return sorted(
        [
            (x, Lipinski.HeavyAtomCount(Chem.MolFromSmiles(x)))
            for x in smiles.split(".")
        ],
        key=lambda x: x[1],
        reverse=True,
    )[0][0]


def sanitize_mol(mol):
    """Sanitized the input molecule

    :param mol: the input molecule
    :return: the sanitized molecule
    """
    s_store_mol = neutralise_charges(
        desalt_compound(Chem.MolToSmiles(mol, isomericSmiles=True))
    )[0]
    store_mol = Chem.MolFromSmiles(s_store_mol)
    if store_mol is None:
        sys.stderr.write(
            "NEUTRALISING MADE NONE MOL"
            + " "
            + s_store_mol
            + " "
            + Chem.MolToSmiles(mol, isomericSmiles=True)
        )
        return None
    return store_mol


def get_path_or_none(xtal_path, xtal, dict_input, dict_key):
    """For each folder with crystal info, and a standard dictionary of expected contents, check whether files exist.

    :param xtal_path: path to crystal folder
    :param xtal: crystal name
    :param dict_input: dictionary containing expected contents of crystal folder
    :param dict_key: file type in dictionary.
    :return: path to file or None if file does not exist.
    """
    if dict_key in dict_input:
        suffix = dict_input[dict_key]
    else:
        print("Key - " + dict_key + " not in dictionary.")
        return None
    path = os.path.join(xtal_path, xtal + suffix)
    if os.path.isfile(path):
        return path
    else:
        print("Path - " + path + " not found.")
        return None


def get_create_target(title):
    """Add a target or return existing target if already exists

    :param title: add a target by title
    :return: the created target
    """
    new_target = Target.objects.get_or_create(title=title)
    print('Target created = ' + str(Target.objects.get_or_create(title=title)[1]))
    return new_target[0]


def add_prot(code, target, xtal_path, xtal, input_dict):
    """Add a protein with a PDB, code and

    :param code: the unique code for this file
    :param target: the target to be linkede to
    :param xtal_path: the path to the crystal directory
    :param xtal: name of the xtal(?)
    :param input_dict: dictionary of files in crystal directory from load_dir
    :return: the created protein
    """

    # Split code by : before the get or create operation and use the first part of the name (split[0])
    # code is normally the xtal directory in the aligned folder, but this may have been modified to have
    # an alternate name added to it - in the form 'directory:alternate_name'.
    code_first_part = code.split(":")[0]
    proteins = Protein.objects.filter(code__contains=code_first_part)
    if proteins.exists():
        new_prot = proteins.first()
    else:
        new_prot = Protein.objects.get_or_create(code=code, target_id=target)
        print('Protein created = ' + str(new_prot[1]))
        new_prot = new_prot[0]

    new_prot.apo_holo = True

    filepaths = {
        'pdb_info': ('pdbs', get_path_or_none(xtal_path, xtal, input_dict, "APO")),
        'bound_info': ('bound', get_path_or_none(xtal_path, xtal, input_dict, "BOUND")),
        'cif_info': ('cifs', get_path_or_none(xtal_path, xtal, input_dict, "CIF")),
        'mtz_info': ('mtzs', get_path_or_none(xtal_path, xtal, input_dict, "MTZ")),
        'map_info': ('maps', get_path_or_none(xtal_path, xtal, input_dict, "PMAP")),
        'sigmaa_info': ('maps', get_path_or_none(xtal_path, xtal, input_dict, "SIGMAA")),
        'diff_info': ('maps', get_path_or_none(xtal_path, xtal, input_dict, "DIFF")),
        'event_info': ('maps', get_path_or_none(xtal_path, xtal, input_dict, "EVENT")),
    }

    to_unpack = {k: v for k, v in filepaths.items() if v[1] is not None}

    for key in to_unpack.keys():
        save_path = os.path.join(to_unpack[key][0], to_unpack[key][1].split('/')[-1])
        path = default_storage.save(save_path, open(to_unpack[key][1], 'rb'))

        setattr(new_prot, key, path)

    new_prot.save()
    return new_prot


def add_projects_to_cmpd(new_comp, projects):
    """Add a project links to a compound

    :param new_comp: the Django compound to add them to
    :param projects:  the list Django projects to add
    :return: the compound with the added projects
    """
    [new_comp.project_id.add(x) for x in projects]
    new_comp.save()
    return new_comp


def calc_cpd(cpd_object, mol, projects):
    """Calculate compound"""

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

    cpd_object.smiles = smiles
    if len(smiles) > Compound._meta.get_field("smiles").max_length:
        print("SMILES TOO LONG")
        return None
    if not len(inchi) > 255:
        cpd_object.inchi = inchi
    else:
        print("INCHI TOO LONG")
        return None
    m = sanitized_mol

    if m is None:
        sys.stderr.write("NONE MOLECULE PRODUCED\n" + smiles + "\n" + inchi)
        return None
    cpd_object.mol_log_p = Chem.Crippen.MolLogP(m)
    cpd_object.mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(m))
    cpd_object.heavy_atom_count = Chem.Lipinski.HeavyAtomCount(m)
    cpd_object.heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(m))
    cpd_object.nhoh_count = Chem.Lipinski.NHOHCount(m)
    cpd_object.no_count = Chem.Lipinski.NOCount(m)
    cpd_object.num_h_acceptors = Chem.Lipinski.NumHAcceptors(m)
    cpd_object.num_h_donors = Chem.Lipinski.NumHDonors(m)
    cpd_object.num_het_atoms = Chem.Lipinski.NumHeteroatoms(m)
    cpd_object.num_rot_bonds = Chem.Lipinski.NumRotatableBonds(m)
    cpd_object.num_val_electrons = Descriptors.NumValenceElectrons(m)
    cpd_object.ring_count = Chem.Lipinski.RingCount(m)
    cpd_object.tpsa = Chem.rdMolDescriptors.CalcTPSA(m)
    # Validate that the compound is unique
    try:
        cpd_object.validate_unique()
        cpd_object.save()
        cpd_object = add_projects_to_cmpd(cpd_object, projects)
        return cpd_object
    except ValidationError:
        if not inchi:
            cpd_object.save()
        cpd_object = Compound.objects.get(inchi=inchi)
        cpd_object = add_projects_to_cmpd(cpd_object, projects)
        return cpd_object


def update_cpd(cpd_id, mol, projects):
    """Update compound"""
    print(mol)
    cpd = cpd_id
    comp = calc_cpd(cpd, mol, projects)
    return comp


def add_comp(mol, projects):
    """Function to add a new compound to the database given an RDKit molecule. Takes an RDKit molecule.

    :param mol: the input RDKit molecule
    :param projects:
    :return: a compound object for the RDKit molecule
    """

    # Now attribute all this meta-deta to the compound object
    new_comp = Compound()
    comp = calc_cpd(new_comp, mol, projects)
    return comp


def add_mol(mol_sd, prot, projects, lig_id="LIG", chaind_id="Z", occupancy=0.0):
    """Function to add a new Molecule to the database

    :param mol_sd: the SDMolBlock of the molecule
    :param prot: the protein it is associated to
    :param projects: the projects it is associated to
    :param lig_id: the 3 letter ligand id
    :param chaind_id: the chain id
    :param occupancy: the occupancy
    :return: the created molecule
    """
    # create mol object from mol_sd
    rd_mol = Chem.MolFromMolFile(mol_sd)

    if rd_mol is None:
        return None

    # See if there is already a molecule with a compound
    old_mols = Molecule.objects.filter(prot_id=prot)

    print('OLD MOLS = ' + str(len(old_mols)))
    # If there's only one
    if len(old_mols) == 1:
        # find the right id (if it exists)
        cpd_id = old_mols[0].cmpd_id
        if cpd_id:
            # update existing compound
            comp_ref = update_cpd(cpd_id, rd_mol, projects)
        else:
            # create new if there's no cpd already
            comp_ref = add_comp(rd_mol, projects)

    else:
        comp_ref = add_comp(rd_mol, projects)

    if comp_ref:
        new_mol = Molecule.objects.get_or_create(prot_id=prot, cmpd_id=comp_ref)[0]
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
    else:
        return None


def parse_proasis(input_string):
    """Parse proasis contact strings

    :param input_string: the Proasis contact string to parse
    :return: a tuple of res_name, res_num, chain_id
    """
    return (
        input_string[:3].strip(),
        int(input_string[5:].strip()),
        input_string[3:5].strip(),
    )


def create_int(prot_res, mol, int_type, interaction):
    """Create a Django interaction object

    :param prot_res: the Django protein residue
    :param mol: the Django molecule
    :param int_type: the interaction type string
    :param interaction: the interaction dictionary
    :return: None
    """
    interation_point = InteractionPoint.objects.get_or_create(
        prot_res_id=prot_res,
        mol_id=mol,
        protein_atom_name=interaction["dstaname"],
        molecule_atom_name=interaction["srcaname"],
    )[0]
    Interaction.objects.get_or_create(
        interaction_version="PR",
        interaction_point=interation_point,
        interaction_type=int_type.get_int_conv("PR", interaction["contactType"]),
        distance=interaction["dis"],
        prot_smarts=interaction["dstType"],
        mol_smarts=interaction["srcType"],
    )


def add_contacts(input_data, target, prot, mol):
    """Add a series of Django contact objects

    :param input_data: the
    :param target: the data - either dict or list - of itneractions
    :param prot: the Django protein object
    :param mol: the Django molecule object
    :return: None
    """
    int_type = IntTypes()
    int_list = []
    if type(input_data) == dict:
        if "results" in input_data:
            int_list = input_data["results"]
    else:
        int_list = input_data
    for interaction in int_list:
        # Ignore Water mediasted hypothesis with Protein for now
        if interaction["hetmoltype"] == "WATER":
            continue
        res_name, res_num, chain_id = parse_proasis(interaction["dstrname"])
        targ_res = TargetResidue.objects.get_or_create(
            target_id=target, res_name=res_name, res_num=res_num, chain_id=chain_id
        )[0]
        prot_res = ProteinResidue.objects.get_or_create(
            targ_res_id=targ_res, prot_id=prot
        )[0]
        create_int(prot_res, mol, int_type, interaction)


def add_map(new_prot, new_target, map_path, map_type):
    """Add a Django map obect

    :param new_prot: the Django protein object
    :param new_target: the Django target object
    :param map_path: the path to the map file
    :param map_type: the two letter code signifyign the type of the map
    :return: the add Django map object
    """
    hotspot_map = HotspotMap.objects.get_or_create(
        map_type=map_type, target_id=new_target, prot_id=new_prot
    )[0]
    hotspot_map.map_info.save(os.path.basename(map_path), File(open(map_path)))
    return hotspot_map


def delete_users(project):
    """Refresh the users for a given project by deleting all of them. Redundant if using iSpyB.

    :param project: the project to remove users from.
    :return: None
    """
    for user_id in project.user_id.all():
        project.user_id.remove(user_id.pk)
    project.save()


def get_create_projects(target, proposal_ref):
    """Add proposals and visits as projects for a given target.

    :param new_target: the target being added
    :param dir_path: the path for where the PROPOSALS and VISITS files are held.
    :return: a list of the projects added.
    """

    # Note that in the loader this is based on information in the PROPOSALS and VISITS files
    # TODO Multiple Visits can be defined in a file apparently - future improvement.
    # TODO NB LIne above in delete_users - redundant if using ISPYB??.
    # For the online loader it comes from the proposal_ref

    projects = []
    # The first word is the ISPY proposal/visit name that is used as the title of the project.
    # It can be set to OPEN in which case there are no users.
    visit = proposal_ref.split()[0]
    project = Project.objects.get_or_create(title=visit)[0]
    projects.append(project)

    # If not open then delete users for the project and re-add them based on supplied fed-ids.
    delete_users(project)

    # Update project_id on target.
    target.project_id.add(project)

    # Remaining words in proposal_ref (if any) must be fedid's which are used to find users information.
    for fedid in proposal_ref.split()[1:]:
        user = User.objects.get_or_create(username=fedid, password="")[0]
        project.user_id.add(user)
    target.save()

    return projects


def remove_not_added(target, xtal_list):
    """Remove any crystals that have not been added this time around. Ensures the database updates, e.g. if someone
    nobody wants a given xtal.

    :param target: the target being checked
    :param xtal_list: a list of folder anmes in the aligned directory
    :return: None
    """
    all_prots = Protein.objects.filter(target_id=target)
    # make sure not to delete any of the computed set proteins (which are protected)
    computed_prots = [mol.pdb for mol in ComputedMolecule.objects.filter(pdb__target_id=target)] 
    unprotected = [x for x in all_prots if x not in computed_prots] 
    
    for prot in unprotected:
        # Code consists of 'directory:alternate_name' if exists (code is renamed based on the metadata)
        code_first_part = prot.code.split(":")[0]
        if code_first_part not in xtal_list:
            prot.delete()
    return None


def save_confidence(mol, file_path, annotation_type="ligand_confidence"):
    """save ligand confidence"""

    input_dict = json.load(open(file_path))
    val_store_dict = ["ligand_confidence_comment", "refinement_outcome", "ligand_confidence_int"]
    for val in val_store_dict:
        if val in input_dict:
            value = input_dict[val]
            if value:
                mol_annot = MolAnnotation.objects.get_or_create(mol_id=mol, annotation_type=annotation_type)[0]
                mol_annot.annotation_text = value
                mol_annot.save()
        else:
            print(val + " not found in " + str(input_dict) + " for mol " + str(mol.prot_id.code))


def load_from_dir(new_target, projects, aligned_path):
    """Load the data for a given target from the directory structure

    param: target (str): the target that has just been created.
    param: project(s) attached to the target
    param: aligned_path (str): the path to the aligned folder for the target.
    return: None

    """

    mols_loaded = 0

    # This seems to be a configuration of the structure of a protein folder
    input_dict = get_dict()

    directories = sorted(os.listdir(aligned_path))
    xtal_list = []

    # Process crystals in the aligned directory
    for xtal in directories:
        mols_loaded += 1
        if not os.path.isdir(os.path.join(aligned_path, xtal)):
            continue
        print(xtal)
        xtal_list.append(xtal)
        xtal_path = os.path.join(aligned_path, xtal)

        pdb_file_path = get_path_or_none(xtal_path, xtal, input_dict, "APO")
        mol_file_path = get_path_or_none(xtal_path, xtal, input_dict, "MOL")
        code = pdb_file_path.split('/')[-1].rsplit("_", 1)[0]

        if pdb_file_path:
            new_prot = add_prot(code=code, target=new_target, xtal_path=xtal_path, xtal=xtal, input_dict=input_dict)
            new_prot.save()
        if mol_file_path:
            new_mol = add_mol(mol_file_path, new_prot, projects)
            if new_mol:
                new_mol.save()

    # Remove proteins for crystals that are not part of the library
    remove_not_added(new_target, xtal_list)

    return mols_loaded


def create_vect_3d(mol, new_vect, vect_ind, vector):
    """Generate the 3D synthesis vectors for a given molecule

    :param mol: the Django molecule object
    :param new_vect: the Django 2d vector object
    :param vect_ind: the index of the vector - since the same 2D vector
    can be different in 3D
    :param vector: the vector coordinates - a 2*3 list of lists.
    :return: None
    """
    if vector:
        new_vect3d = Vector3D.objects.get_or_create(
            mol_id=mol, vector_id=new_vect, number=vect_ind
        )[0]
        # The start position
        new_vect3d.start_x = float(vector[0][0])
        new_vect3d.start_y = float(vector[0][1])
        new_vect3d.start_z = float(vector[0][2])
        # The end position
        new_vect3d.end_x = float(vector[1][0])
        new_vect3d.end_y = float(vector[1][1])
        new_vect3d.end_z = float(vector[1][2])
        new_vect3d.save()


def get_vectors(mols):
    """Get the vectors for a given molecule

    :param mols: the Django molecules to get them from
    :return: None
    """
    vect_types = VectTypes()
    for mol in mols:
        if "." in mol.smiles:
            print("SKIPPING - FRAGMENT: " + str(mol.smiles))
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
                new_vect = Vector.objects.get_or_create(
                    smiles=smiles, cmpd_id=mol.cmpd_id, type=vect_choice
                )[0]
                create_vect_3d(mol, new_vect, vect_ind, vectors[vect_type][vector])


def cluster_mols(rd_mols, mols, target):
    """Cluster a series of 3D molecules

    :param rd_mols: the RDKit molecules to cluster
    :param mols:  the Django moleculs they refer to
    :param target:  the Django target it refers to
    :return: None
    """
    id_mols = [x.pk for x in mols]
    out_data = run_lig_cluster(rd_mols, id_mols)
    for clust_type in out_data:
        for cluster in out_data[clust_type]:
            # look for molgroup with same coords - need to implement tolerance?
            mol_group = search_for_molgroup_by_coords(coords=[out_data[clust_type][cluster]["centre_of_mass"][0],
                                                              out_data[clust_type][cluster]["centre_of_mass"][1],
                                                              out_data[clust_type][cluster]["centre_of_mass"][2]],
                                                      target=target.title)
            if not mol_group:
                mol_group = MolGroup()
            if clust_type != "c_of_m":
                mol_group.group_type = "PC"
            else:
                mol_group.group_type = "MC"
            mol_group.target_id = target
            mol_group.x_com = out_data[clust_type][cluster]["centre_of_mass"][0]
            mol_group.y_com = out_data[clust_type][cluster]["centre_of_mass"][1]
            mol_group.z_com = out_data[clust_type][cluster]["centre_of_mass"][2]
            mol_group.description = clust_type
            mol_group.save()
            for mol_id in out_data[clust_type][cluster]["mol_ids"]:
                if mol_id not in [a['id'] for a in mol_group.mol_id.values()]:
                    this_mol = Molecule.objects.get(id=mol_id)
                    mol_group.mol_id.add(this_mol)


def centre_of_points(list_of_points):
    """average list of points"""

    cp = np.average(list_of_points, axis=0)
    return cp


def centre_of_mass(mol):
    """calculate centre of mass"""

    numatoms = mol.GetNumAtoms()
    conf = mol.GetConformer()
    if not conf.Is3D():
        return 0
    # get coordinate of each atoms
    pts = np.array([list(conf.GetAtomPosition(atmidx)) for atmidx in range(numatoms)])
    atoms = [atom for atom in mol.GetAtoms()]
    mass = Descriptors.MolWt(mol)
    # get center of mass
    center_of_mass = np.array(np.sum(atoms[i].GetMass() * pts[i] for i in range(numatoms))) / mass
    return center_of_mass


def process_site(rd_mols):
    """process site"""

    coms = [centre_of_mass(mol) for mol in rd_mols]
    centre = centre_of_points(coms)
    print('CENTRE: ' + str(centre))
    return centre


def get_coord_limits(coord):
    """get cooordinate limits"""

    lower_limit = float('.'.join([str(coord).split('.')[0], str(coord).split('.')[1][:2]]))
    if lower_limit > 0:
        upper_limit = lower_limit + 0.01
    else:
        tmp = lower_limit - 0.01
        upper_limit = lower_limit
        lower_limit = tmp
    return lower_limit, upper_limit


def search_for_molgroup_by_coords(coords, target):
    """search for a molgroup by list of coordinates"""

    x = coords[0]
    y = coords[1]
    z = coords[2]

    limit_list = []

    for coord in x, y, z:
        lower, upper = get_coord_limits(coord)
        limit_list.append([lower, upper])

    search = MolGroup.objects.filter(target_id__title=target, x_com__gte=limit_list[0][0], x_com__lte=limit_list[0][1],
                                     y_com__gte=limit_list[1][0], y_com__lte=limit_list[1][1],
                                     z_com__gte=limit_list[2][0],
                                     z_com__lte=limit_list[2][1])

    if len(search) == 1:
        mol_group = search[0]
    else:
        return None

    return mol_group


def search_for_molgroup_by_description(description, target):
    """search for a molgroup by description"""

    search = MolGroup.objects.filter(target_id__title=target, description=description)
    print(str('matching_sites = ')+str(len(search)))
    if len(search) == 1:
        mol_group = search[0]

    elif len(search) > 1:
        for molgroup in search:
            molgroup.delete()
        return None
    else:
        return None

    return mol_group


def analyse_mols(mols, target, specified_site=False, site_description=None):
    """Analyse a list of molecules for a given target

    :param mols: the Django molecules to analyse
    :param target: the Django target
    :param specified_site:
    :param site_description:
    :return: None
    """
    rd_mols = [Chem.MolFromMolBlock(x.sdf_info) for x in mols]
    if not specified_site:

        cluster_mols(rd_mols, mols, target)

    else:

        centre = process_site(rd_mols)

        # look for molgroup with same target and description
        mol_group = search_for_molgroup_by_description(target=target.title, description=site_description)

        if not mol_group:
            mol_group = MolGroup()
        mol_group.group_type = "MC"
        mol_group.target_id = target
        mol_group.x_com = centre[0]
        mol_group.y_com = centre[1]
        mol_group.z_com = centre[2]
        mol_group.description = site_description
        mol_group.save()

        ids = [m.id for m in mols]

        print([a['id'] for a in mol_group.mol_id.values()])

        for mol_id in ids:
            if mol_id not in [a['id'] for a in mol_group.mol_id.values()]:
                print(mol_id)
                this_mol = Molecule.objects.get(id=mol_id)
                mol_group.mol_id.add(this_mol)

    get_vectors(mols)


def rename_proteins(names_csv):
    """rename proteins based on the contents of alternatenames.csv (which has been created from metdata.csv """

    names_frame = pd.read_csv(names_csv)

    for _, row in names_frame.iterrows():
        mol_target = row['name']
        alternate_name = row['alternate_name']
        # Remove the replacement of '_0' - this was inconsistently applied as some folders are '_1'
        # The Protein code will be modified to be of format 'xtal_directory:alternate_name'
        new_name = str(mol_target).strip() + ':' + str(alternate_name).strip()

        prots = Protein.objects.filter(code=mol_target)
        for prot in prots:
            print('changing prot name to: ' + new_name)
            prot.code = new_name
            prot.save()


def relative_to_media_root(filepath, media_root=settings.MEDIA_ROOT):
    """Calculate a relative path from the file path to the media root"""
    relative_path = os.path.relpath(filepath, media_root)
    return relative_path


def analyse_target(target_name, aligned_path):
    """Analyse all the molecules for a particular target.

    param: target_name (str): the string title of the target. This will uniquely identify it.
    param: aligned_path (str): the path to the data in the aligned directory for the target.
    return: None

    """

    # Get Target from database
    target = Target.objects.get(title=target_name)
    target.root_data_directory = relative_to_media_root(aligned_path)
    target.save()

    mols = list(Molecule.objects.filter(prot_id__target_id=target))

    # This can probably be improved to count molecules as they are processed when the code is further refactored
    mols_processed = len(mols)

    print("Analysing " + str(len(mols)) + " molecules for " + target_name)

    # Do site mapping
    if os.path.isfile(os.path.join(aligned_path, 'metadata.csv')):

        target.metadata.save(
            os.path.basename(os.path.join(aligned_path, 'metadata.csv')),
            File(open(os.path.join(aligned_path, 'metadata.csv')))
        )

        # remove any existing files so that we don't create a messy file when appending
        if os.path.isfile(os.path.join(aligned_path, 'hits_ids.csv')):
            os.remove(os.path.join(aligned_path, 'hits_ids.csv'))

        if os.path.isfile(os.path.join(aligned_path, 'sites.csv')):
            os.remove(os.path.join(aligned_path, 'sites.csv'))

        if os.path.isfile(os.path.join(aligned_path, 'alternate_names.csv')):
            os.remove(os.path.join(aligned_path, 'alternate_names.csv'))

        new_frame = pd.read_csv(os.path.join(aligned_path, 'metadata.csv'))
        new_frame.sort_values(by='site_name', inplace=True)

        # one file for new names
        with open(os.path.join(aligned_path, 'alternate_names.csv'), 'a') as f:
            f.write('name,alternate_name\n')

            for _, row in new_frame.iterrows():
                if isinstance(row['alternate_name'], str):
                    crystal_name = row['crystal_name']
                    # find the correct crystal
                    crystal = Protein.objects.filter(code__contains=crystal_name, target_id=target)
                    alternate_name = row['alternate_name']
                    # Only take first part of code
                    for crys in list(set([c.code.split(":")[0] for c in crystal])):
                        f.write(str(crys) + ',' + str(alternate_name) + '\n')

        # hits and sites files
        site_mapping = {}
        unique_sites = list(set(list(new_frame['site_name'])))
        for i in range(0, len(sorted(unique_sites))):
            site_mapping[unique_sites[i]] = i

        with open(os.path.join(aligned_path, 'hits_ids.csv'), 'a') as f:
            f.write('crystal_id,site_number\n')

            for _, row in new_frame.iterrows():
                crystal_name = row['crystal_name']
                crystal = Protein.objects.filter(code__contains=crystal_name, target_id=target)
                site = row['site_name']
                s_id = site_mapping[site]
                for crys in list(set([c.code for c in crystal])):
                    f.write(str(crys) + ',' + str(s_id) + '\n')

        with open(os.path.join(aligned_path, 'sites.csv'), 'a') as f:
            f.write('site,id\n')
            for key in site_mapping.keys():
                f.write(str(key) + ',' + str(site_mapping[key]) + '\n')

    if os.path.isfile(os.path.join(aligned_path, 'hits_ids.csv')) and os.path.isfile(
            os.path.join(aligned_path, 'sites.csv')):

        hits_sites = pd.read_csv(os.path.join(aligned_path, 'hits_ids.csv'))
        sites = pd.read_csv(os.path.join(aligned_path, 'sites.csv'))
        sites.sort_values(by='site', inplace=True)

        # delete the old molgroups first
        mgs = MolGroup.objects.filter(target_id=target)
        for m in mgs:
            m.delete()

        for _, row in sites.iterrows():
            description = row['site']
            number = row['id']
            print('Processing user input site: ' + str(description))
            matches = []
            for _, row in hits_sites.iterrows():
                if str(row['site_number']) == str(number):
                    matches.append(row['crystal_id'])
            print('HIT IDS: ' + str(matches))
            print('\n')
            if matches:
                mols = list(Molecule.objects.filter(prot_id__target_id=target, prot_id__code__in=matches))
                analyse_mols(mols=mols, target=target, specified_site=True, site_description=description)

    if os.path.isfile(os.path.join(aligned_path, 'alternate_names.csv')):
        rename_proteins(names_csv=os.path.join(aligned_path, 'alternate_names.csv'))
    else:
        analyse_mols(mols=mols, target=target)

    # move anything that's not a directory in 'aligned' up a level
    files = (f for f in os.listdir(aligned_path)
             if os.path.isfile(os.path.join(aligned_path, f)))

    for f in files:
        shutil.move(os.path.join(aligned_path, f), os.path.join(aligned_path, f).replace('aligned', ''))

    # delete NEW_DATA VISITS PROPOSALS. These are not used by the new loader but might be in old data sets.
    to_delete = ['NEW_DATA', 'VISITS', 'PROPOSALS']
    for file in to_delete:
        filepath = os.path.join(aligned_path.replace('aligned', ''), file)
        if os.path.isfile(filepath):
            os.remove(filepath)

    # last step - zip up the input file and move it to the archive
    zipped = shutil.make_archive(aligned_path.replace('aligned', ''), 'zip', aligned_path.replace('aligned', ''))
    # shutil.move(zipped, os.path.join(settings.MEDIA_ROOT, 'targets', os.path.basename(zipped)))
    target.zip_archive.name = relative_to_media_root(zipped)
    target.save()

    return mols_processed


def process_target(new_data_folder, target_name, proposal_ref):
    """Process the full target dataset.

    :param target new_data_folder: path where the extracted target folder is located
    :param target_name: Name of the target - should be the same as folder with the tmp_folder
    :param proposal_ref: A reference to the proposal/visit used for connecting the target to a project/users
    :return: mols_loaded, mols_processed
    """
    mols_loaded = 0
    mols_processed = 0

    # e.g. code/media/tmp/new_data -> replaced by the target within the directory: code/media/tmp/new_data/Mpro .
    target_path = os.path.join(new_data_folder, target_name)

    # path to save the media to
    # /code/media/
    media_root = settings.MEDIA_ROOT
    # /code/media/targets
    upload_path = os.path.join(media_root, 'targets')
    os.makedirs(upload_path, exist_ok=True)

    # e.g. /code/media/targets/mArh
    target_upload_path = os.path.join(upload_path, target_name)

    # if target data already exists, remove existing data in the media directory for the target
    if os.path.isdir(target_upload_path):
        shutil.rmtree(target_upload_path)

    print('Saving uploaded data to ' + upload_path)
    # move the whole folder from the upload directory to the media directory
    # This creates the initial data in the aligned directory.
    shutil.move(target_path, upload_path)

    # change the target_path to the new 'aligned' directory
    aligned_path = os.path.join(upload_path, target_name, 'aligned')

    print('ALIGNED_PATH: ' + aligned_path)
    # Check if there is any data to process
    if os.path.isdir(aligned_path):
        # Create the target if required
        new_target = get_create_target(target_name)

        # Create a project attached to the target with proposal/visit information if it exists.
        projects = get_create_projects(new_target, proposal_ref)

        # "Upsert" molecule data from the aligned folder for the target
        mols_loaded = load_from_dir(new_target, projects, aligned_path)

        # This updates files like metadata.csv, alternatename.csv, hits_ids.csv etc.
        if mols_loaded:
            mols_processed = analyse_target(target_name, aligned_path)
    else:
        print("Aligned folder is missing - no data to add: " + aligned_path)

    return mols_loaded, mols_processed


def add_tset_warning(validate_dict, location, error, line_number):

    validate_dict['Location'].append(location)
    validate_dict['Error'].append(error)
    validate_dict['Line number'].append(line_number)
    return validate_dict

def check_meatadata_row(validated, input_validate_dict, row, idx):
    """Validate the metadata.csv file to check basic formatting is correct

    Metedata.csv has the following columns:
    crystal_name: must not be spaces or null and should contain the RealCrystalName
    RealCrystalName: must not be spaces or null
    smiles: must not be null
    new_smiles: no specific validation
    alternate_name: no specific validation
    site_name: whole column should either be null or not null (no partial columns)
    pdb_entry: no specific validation

    :return:
        validated: boolean (updated)
        validate_dict: (updated)
    """

    if row['RealCrystalName'].isspace() or row['RealCrystalName'] == 'nan':
        add_tset_warning(input_validate_dict, 'Metadata.csv', 'RealCrystalName spaces or null', idx + 2)
        validated = False
    if row['crystal_name'].isspace() or row['RealCrystalName'] == 'nan':
        add_tset_warning(input_validate_dict, 'Metadata.csv', 'Crystal name spaces or null', idx + 2)
        validated = False
    if row['RealCrystalName'] not in row['crystal_name']:
        add_tset_warning(input_validate_dict, 'Metadata.csv', 'Crystal name does not contain RealCrystalName', idx + 2)
        validated = False
    if row['smiles'] == 'nan':
        add_tset_warning(input_validate_dict, 'Metadata.csv', 'Smiles null', idx + 2)
        validated = False

    return validated, input_validate_dict


def check_metadata(metadata_file, input_validate_dict):
    """Validate the metadata.csv file to check basic formatting is correct

    :return:
        validated: boolean
        validate_dict:
    """
    validated = True
    # Metedata.csv has the following columns:
    # crystal_name: must not be spaces or null and should contain the RealCrystalName
    # RealCrystalName: must not be spaces or null
    # smiles: must not be null
    # new_smiles: no specific validation
    # alternate_name: no specific validation
    # site_name: whole column should either be null or not null (no partial columns)
    # pdb_entry: no specific validation

    meta_dataframe = pd.read_csv(metadata_file)

    # File level checks.
    meta_sites = meta_dataframe['site_name']
    if meta_sites.isnull().values.all() or meta_sites.notnull().values.all():
        pass
    else:
        add_tset_warning(input_validate_dict, 'Metadata.csv',
                    'site_name column should either be completely filled or completely null', 0)
        validated = False

    meta_dataframe['crystal_name'] = meta_dataframe['crystal_name'].astype(str)
    meta_dataframe['RealCrystalName'] = meta_dataframe['RealCrystalName'].astype(str)
    meta_dataframe['smiles'] = meta_dataframe['smiles'].astype(str)

    # Loop through metadata doing basic checks on each row
    for idx, (_, row) in enumerate(meta_dataframe.iterrows()):
        validated, input_validate_dict = check_meatadata_row(validated, input_validate_dict, row, idx)

    return validated, input_validate_dict


def validate_target(new_data_folder, target_name, proposal_ref):
    """Validate the target dataset. This will initially just be structurally

    :param target new_data_folder: path where the extracted target folder is located
    :param target_name: Name of the target - should be the same as folder with the tmp_folder
    :param proposal_ref: A reference to the proposal/visit used for connecting the target to a project/users
    :return:
    """
    validate_dict = {'Location': [], 'Error': [], 'Line number': []}

    # Check if there is any data to process
    target_path = os.path.join(new_data_folder, target_name)

    if not os.path.isdir(target_path):
        validate_dict = add_tset_warning(validate_dict, 'Folder', f'No folder matching target name in extracted zip file {new_data_folder}, {target_name}', 0)

    aligned_path = os.path.join(target_path, 'aligned')

    if not os.path.isdir(aligned_path):
        validate_dict = add_tset_warning(validate_dict, 'Folder', 'No aligned folder present in target name folder', 0)

    # Check if there is a metadata.csv file to process
    metadata_file = os.path.join(aligned_path, 'metadata.csv')
    if os.path.isfile(metadata_file):
        validated, validate_dict = check_metadata(metadata_file, validate_dict)
    else:
        validate_dict = add_tset_warning(validate_dict, 'File', 'No metedata.csv file present in the aligned folder', 0)
        validated = False

    return validated, validate_dict
