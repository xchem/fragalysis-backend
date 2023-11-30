"""
Target configuration functions for uploading a Target Dataset to Fragalysis.
This was originally sourced from loader/config.py
"""


def _mol_choices():
    """
    Get the data for the molecule choices
    :return:
    """
    PROASIS = "PR"
    SDF = "SD"
    HYDROGEN = "HA"
    HYDROGEN_AM1_CHARGES = "HC"
    mol_choices = (
        (PROASIS, "Proasis molecule", ".mol", "MOL"),
        (SDF, "Sdf molecule", ".sdf", "SDF"),
        (HYDROGEN, "Hydrogens added ", "_h.mol", "H_MOL"),
        (
            HYDROGEN_AM1_CHARGES,
            "Mol2 format with Hydrogens and AM1 BCC",
            ".mol2",
            "MOL2",
        ),
    )
    return mol_choices, PROASIS, SDF


def _prot_choices():
    """
    Get the data for the protein choices
    :return:
    """
    APO = "AP"
    STRIPPED = "ST"
    TLEAPED = "TL"
    CHUNKED = "CH"
    BOUND = "BO"
    prot_choices = (
        (APO, "Apo", "_apo.pdb", "APO"),
        (STRIPPED, "Stripped", "_no_buffer_altlocs.pdb", "STRIPPED"),
        (TLEAPED, "Tleaped", "_tleap.pdb", "TLEAP"),
        (CHUNKED, "Chunked", "_chunk.pdb", "CHUNK"),
        (BOUND, "Bound", '_bound.pdb', "BOUND"),
    )
    return prot_choices, APO


def _djangofy_choices(inputs):
    """
    Convert a list of tuples into the required format
    :param input:
    :return:
    """
    return [x[:2] for x in inputs]


def get_dict():
    """
    Convert the mol and protein choices into a dictionary
    :return:
    """
    out_d = {}
    for prot in _prot_choices()[0]:
        out_d[prot[3]] = prot[2]
    for mol in _mol_choices()[0]:
        out_d[mol[3]] = mol[2]
    out_d["BOUND"] = "_bound.pdb"
    out_d["EVENT"] = "_event.ccp4"
    out_d["SIGMAA"] = "_2fofc.map"
    out_d["DIFF"] = "_fofc.map"
    # TODO: Look up how this works and allow either ccp4 or map
    # out_d["EVENT"] = "_event.ccp4"
    # out_d["SIGMAA"] = "_2fofc.ccp4"
    # out_d["DIFF"] = "_fofc.ccp4"
    out_d["MTZ"] = ".mtz"
    # out_d["PMAP"] = "_event.ccp4"
    out_d["PPDB"] = "_pandda.pdb"
    out_d["PJSON"] = "_pandda.json"
    out_d["PMTZ"] = "_pandda.mtz"
    # optional ones - contacts and hotspots
    out_d["CONTACTS"] = "_contacts.json"
    out_d["CONFIDENCE"] = "_lig_conf.json"
    out_d["ACC"] = "_acceptor.ccp4"
    out_d["DON"] = "_donor.ccp4"
    out_d["LIP"] = "_apolar.ccp4"
    out_d["TRANS"] = "_transform.json"
    out_d["DESOLV"] = "_apo-desolv.pdb"
    out_d["HEADER"] = "_header.pdb"
    return out_d


def get_mol_choices():
    """
    Get the molecule choices for django
    :return:
    """
    return _djangofy_choices(_mol_choices()[0]), _mol_choices()[1]


def get_prot_choices():
    """
    Get the protein choices for django
    :return:
    """
    return _djangofy_choices(_prot_choices()[0]), _prot_choices()[1]
