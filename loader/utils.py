from loader.functions import NeutraliseCharges,desalt_compound

def sanitize_mol(mol):
    """
    Sanitized the input molecule
    :param mol: the input molecule
    :return: the sanitized molecule
    """
    s_store_mol = NeutraliseCharges(desalt_compound(Chem.MolToSmiles(mol, isomericSmiles=True)))[0]
    store_mol = Chem.MolFromSmiles(s_store_mol)
    if store_mol is None:
        sys.stderr.write("NEUTRALISING MADE NONE MOL" + " " + s_store_mol + " " + Chem.MolToSmiles(mol, isomericSmiles=True))
        return None
    return store_mol

def get_path_or_none(new_path,xtal,suffix):
    """
    Get a path or none - for loader
    :param new_path:
    :param xtal:
    :param suffix:
    :return:
    """
    path = os.path.join(new_path, xtal + suffix)
    if os.path.isfile(path):
        return path
    else:
        return None