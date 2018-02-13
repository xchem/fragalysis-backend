
from rdkit import Chem
from rdkit.Chem import Lipinski,Descriptors,AllChem

# Contribution to the RDKit from Hans de Winter
def _InitialiseNeutralisationReactions():
    """Contribution from Hans de Winter"""
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None


def NeutraliseCharges(smiles, reactions=None):
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
        return (Chem.MolToSmiles(mol, True), True)
    else:
        return (smiles, False)

def desalt_compound(smiles):
    """Function to desalt compound a given smiles string
    Takes a smiles string
    Returns a desalted smiles string."""
    # Chose the biggest fragment, after splitting into fragments
    return sorted([(x, Lipinski.HeavyAtomCount(Chem.MolFromSmiles(x))) for x in smiles.split(".")], key=lambda x: x[1], reverse=True)[0][0]


# def generate_data_for_smis(smi_list):
#     mols = [Chem.MolFromSmiles(x) for x in smi_list]
#     mols = [x for x in mols if x is not None]
#     for m in mols:
#         # Now generate a series of data for these molecules