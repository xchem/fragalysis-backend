"""
target_setupload.py Functions for uploading a Target Dataset to Fragalysis.

This was originally sourced from the code in the fragment-loader repo:
loaders.py
functions.py
"""
import logging

# import json
import os
import sys
from typing import Any, Dict

# from viewer.target_set_config import get_dict
import numpy as np
import pandas as pd
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.core.files import File
from frag.alysis.run_clustering import run_lig_cluster
from frag.network.decorate import get_3d_vects_for_mol
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski

from hotspots.models import HotspotMap
from hypothesis.definitions import IntTypes, VectTypes
from hypothesis.models import (
    Interaction,
    InteractionPoint,
    TargetResidue,
    Vector,
    Vector3D,
)
from scoring.models import SiteObservationGroup
from viewer.models import (
    Compound,
    Project,
    SiteObservation,
    SiteObservationTag,
    TagCategory,
    Target,
)

# import shutil
# import datetime


logger = logging.getLogger(__name__)

_reactions = None


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


def neutralise_charges(smiles, reactions=None):
    """Contribution from Hans de Winter"""
    global _reactions  # pylint: disable=global-statement
    if reactions is None:
        if _reactions is None:
            _reactions = _InitialiseNeutralisationReactions()
        reactions = _reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for _, (reactant, product) in enumerate(reactions):
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
    :return: the molecule
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
        logger.warning("Key - '%s' not in dictionary", dict_key)
        return None
    path = os.path.join(xtal_path, xtal + suffix)
    if os.path.isfile(path):
        return path
    else:
        logger.warning("Path - '%s' not found", path)
        return None


def get_create_target(title):
    """Add a target or return existing target if already exists

    :param title: add a target by title
    :return: the created target
    """
    new_target = Target.objects.get_or_create(title=title)
    logger.debug(
        "Target created new_target='%s'", Target.objects.get_or_create(title=title)[1]
    )
    return new_target[0]


def add_projects_to_cmpd(new_comp, projects):
    """Add a project links to a compound

    :param new_comp: the Django compound to add them to
    :param projects:  the list Django projects to add
    :return: the compound with the added projects
    """
    for project in projects:
        new_comp.project_id.add(project)
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
    len_smiles = len(smiles)
    # Rather than using the protected member,
    # why not introduce constants in the model?
    if (
        len_smiles
        > Compound._meta.get_field(  # pylint: disable=protected-access
            "smiles"
        ).max_length
    ):
        logger.warning("SMILES too long (%s) [%d]", smiles, len_smiles)
        return None
    len_inchi = len(inchi)
    if len_inchi <= 255:
        cpd_object.inchi = inchi
    else:
        logger.warning("INCHI too long (%s) [%d]", inchi, len_inchi)
        return None

    m = sanitized_mol
    if m is None:
        msg = "NONE MOLECULE PRODUCED\n" + smiles + "\n" + inchi
        sys.stderr.write(msg)
        logger.warning(msg)
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
    logger.debug("update_cpd(%s, %s, ...)", cpd_id, mol)
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


def create_int(target_res, site_obvs, int_type, interaction):
    """Create a Django interaction object

    :param target_res: the Django target protein residue
    :param site_obvs: the Django site_observation
    :param int_type: the interaction type string
    :param interaction: the interaction dictionary
    :return: None
    """
    interation_point = InteractionPoint.objects.get_or_create(
        target_res=target_res,
        site_observation=site_obvs,
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


def add_contacts(input_data, target, site_obvs):
    """Add a series of Django contact objects

    :param input_data: the
    :param target: the data - either dict or list - of itneractions
    :param site_obvs: the Django site_observation object
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
        create_int(targ_res, site_obvs, int_type, interaction)


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
    hotspot_map.map_info.save(
        os.path.basename(map_path), File(open(map_path, encoding='utf-8'))
    )
    return hotspot_map


def delete_users(project):
    """Refresh the users for a given project by deleting all of them. Redundant if using iSpyB.

    :param project: the project to remove users from.
    :return: None
    """
    for user_id in project.user_id.all():
        project.user_id.remove(user_id.pk)
    project.save()


def get_create_projects(target, proposal_ref, proposal_code='lb'):
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
    # The first word is the ISPyB proposal/visit name.
    # This is used as the title of the project (e.g. "lb12345-4")
    visit = proposal_ref.split()[0]
    # If the visit is not prefixed by the proposal code
    # (typically the 2-letter sequence "lb") then prefix it.
    if visit[0].isdigit():
        visit = f"{proposal_code}{visit}"
    project = Project.objects.get_or_create(title=visit)[0]
    projects.append(project)

    # If not open then delete users for the project and re-add them based on supplied fed-ids.
    delete_users(project)

    # Update project on target.
    target.project = project

    # The remaining words in proposal_ref (if any)
    # are expected to be fedid's (user IDs) which are used to find user information.
    num_users = 0
    for fedid in proposal_ref.split()[1:]:
        user = User.objects.get_or_create(username=fedid, password="")[0]
        project.user_id.add(user)
        num_users += 1
    if num_users == 0:
        project.open_to_public = True

    target.upload_progess = 10.00
    target.save()

    return projects


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
            logger.debug("SKIPPING - FRAGMENT: %s", mol.smiles)
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


def search_for_molgroup_by_coords(coords, target):
    """search for a molgroup by list of coordinates"""

    x = coords[0]
    y = coords[1]
    z = coords[2]

    limit_list = []

    for coord in x, y, z:
        lower, upper = get_coord_limits(coord)
        limit_list.append([lower, upper])

    search = SiteObservationGroup.objects.filter(
        target_id__title=target,
        x_com__gte=limit_list[0][0],
        x_com__lte=limit_list[0][1],
        y_com__gte=limit_list[1][0],
        y_com__lte=limit_list[1][1],
        z_com__gte=limit_list[2][0],
        z_com__lte=limit_list[2][1],
    )

    if len(search) == 1:
        mol_group = search[0]
    else:
        return None

    return mol_group


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
            mol_group = search_for_molgroup_by_coords(
                coords=[
                    out_data[clust_type][cluster]["centre_of_mass"][0],
                    out_data[clust_type][cluster]["centre_of_mass"][1],
                    out_data[clust_type][cluster]["centre_of_mass"][2],
                ],
                target=target.title,
            )
            if not mol_group:
                mol_group = SiteObservationGroup()
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
                    this_mol = SiteObservation.objects.get(id=mol_id)
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
    center_of_mass = (
        np.array(np.sum(atoms[i].GetMass() * pts[i] for i in range(numatoms))) / mass
    )
    return center_of_mass


def calc_site_centre(rd_mols):
    """Calculate the centre of the site's molecules based on the centre of mass"""

    coms = [centre_of_mass(mol) for mol in rd_mols]
    centre = centre_of_points(coms)
    logger.debug('CENTRE: %s', centre)
    return centre


def get_coord_limits(coord):
    """get cooordinate limits"""

    lower_limit = float(
        '.'.join([str(coord).split('.')[0], str(coord).split('.')[1][:2]])
    )
    if lower_limit > 0:
        upper_limit = lower_limit + 0.01
    else:
        tmp = lower_limit - 0.01
        upper_limit = lower_limit
        lower_limit = tmp
    return lower_limit, upper_limit


def search_for_site_obvsgroup_by_description(description, target):
    """search for a molgroup by description"""

    search = SiteObservationGroup.objects.filter(
        target_id__title=target, description=description
    )
    logger.debug("len(search)=%d", search.count())
    if search.count() == 1:
        site_obvs_group = search[0]

    elif search.count() > 1:
        search.delete()
        return None
    else:
        return None

    return site_obvs_group


def specifc_site(rd_mols, site_observations, target, site_description=None):
    """Update/Create mol_groups and molecule_tags with site information
    :param rd_mols: the molecules to add to the site (rd form)
    :param site_observations: the site_observations to add to the site
    :param target: the Django target
    :param site_description:
    :return: None
    """

    # look for molgroup with same target and description
    site_obvs_group = search_for_site_obvsgroup_by_description(
        target=target.title, description=site_description
    )

    if not site_obvs_group:
        site_obvs_group = SiteObservationGroup()

    site_obvs_group.group_type = "MC"
    site_obvs_group.target_id = target
    centre = calc_site_centre(rd_mols)
    site_obvs_group.x_com = centre[0]
    site_obvs_group.y_com = centre[1]
    site_obvs_group.z_com = centre[2]
    site_obvs_group.description = site_description
    site_obvs_group.save()

    # A molecule tag record may exist already, but won't the first time the
    # target is loaded.

    try:
        site_obvs_tag = SiteObservationTag.objects.get(
            upload_name=site_description, target_id=target.id
        )
    except SiteObservationTag.DoesNotExist:
        site_obvs_tag = None

    if not site_obvs_tag:
        # New site/tag or the tag has been deleted
        site_obvs_tag = SiteObservationTag()
        site_obvs_tag.tag = site_description
        site_obvs_tag.upload_name = site_description
        site_obvs_tag.category = TagCategory.objects.get(category='Sites')
        site_obvs_tag.target = target
        site_obvs_tag.mol_group = site_obvs_group
        site_obvs_tag.save()
    else:
        # Tag already exists
        # Apart from the new mol_group and molecules, we shouldn't be
        # changing anything.
        site_obvs_tag.mol_group = site_obvs_group
        site_obvs_tag.save()

    for site_obvs in site_observations:
        if site_obvs not in site_obvs_group.site_observation:
            logger.debug("site_obvs_group site_obvs_id=%s", site_obvs.id)
            site_obvs_group.site_observation.add(site_obvs)

        if site_obvs not in site_obvs_tag.site_observations:
            logger.debug("site_obvs_tag site_obvs_id=%s", site_obvs.id)
            site_obvs_tag.site_observations.add(site_obvs)


def analyse_mols(mols, target, specified_site=False, site_description=None):
    """Check if molecules belong to a cluster or a specific site for a given target

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
        specifc_site(rd_mols, mols, target, site_description)

    get_vectors(mols)


def relative_to_media_root(filepath, media_root=settings.MEDIA_ROOT):
    """Calculate a relative path from the file path to the media root"""
    relative_path = os.path.relpath(filepath, media_root)
    return relative_path


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
        add_tset_warning(
            input_validate_dict,
            'Metadata.csv',
            'RealCrystalName spaces or null',
            idx + 2,
        )
        validated = False
    if row['crystal_name'].isspace() or row['RealCrystalName'] == 'nan':
        add_tset_warning(
            input_validate_dict, 'Metadata.csv', 'Crystal name spaces or null', idx + 2
        )
        validated = False
    if row['RealCrystalName'] not in row['crystal_name']:
        add_tset_warning(
            input_validate_dict,
            'Metadata.csv',
            'Crystal name does not contain RealCrystalName',
            idx + 2,
        )
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
        add_tset_warning(
            input_validate_dict,
            'Metadata.csv',
            'site_name column should either be completely filled or completely null',
            0,
        )
        validated = False

    meta_dataframe['crystal_name'] = meta_dataframe['crystal_name'].astype(str)
    meta_dataframe['RealCrystalName'] = meta_dataframe['RealCrystalName'].astype(str)
    meta_dataframe['smiles'] = meta_dataframe['smiles'].astype(str)

    # Loop through metadata doing basic checks on each row
    for idx, (_, row) in enumerate(meta_dataframe.iterrows()):
        validated, input_validate_dict = check_meatadata_row(
            validated, input_validate_dict, row, idx
        )

    return validated, input_validate_dict


def validate_target(new_data_folder, target_name, proposal_ref):
    """Validate the target dataset. This will initially just be structurally

    :param target new_data_folder: path where the extracted target folder is located
    :param target_name: Name of the target - should be the same as folder with the tmp_folder
    :param proposal_ref: A reference to the proposal/visit used for connecting the target to a project/users
    :return:
    """
    # Don't need
    del proposal_ref

    validate_dict: Dict[str, Any] = {'Location': [], 'Error': [], 'Line number': []}

    # Check if there is any data to process
    target_path = os.path.join(new_data_folder, target_name)

    # Assume success...
    validated = True

    # A target directory must exist
    if not os.path.isdir(target_path):
        validate_dict = add_tset_warning(
            validate_dict,
            'Folder',
            'Folder does not match target name.'
            f' Expected "{target_name}".'
            f' Is the upload called "{target_name}.zip"?',
            0,
        )
        # No point in checking anything else if this check fails
        validated = False

    if validated:
        # An 'aligned' directory must exist
        aligned_path = os.path.join(target_path, 'aligned')
        if not os.path.isdir(aligned_path):
            validate_dict = add_tset_warning(
                validate_dict,
                'Folder',
                'No aligned folder present.'
                f' Expected "{target_name}/{aligned_path}"',
                0,
            )

    if validated:
        # A metadata.csv file must exist
        metadata_file = os.path.join(aligned_path, 'metadata.csv')
        if os.path.isfile(metadata_file):
            validated, validate_dict = check_metadata(metadata_file, validate_dict)
        else:
            validate_dict = add_tset_warning(
                validate_dict,
                'File',
                'No metedata file present.'
                f' Expected "{target_name}/{aligned_path}/{metadata_file}"',
                0,
            )
            validated = False

    return validated, validate_dict
