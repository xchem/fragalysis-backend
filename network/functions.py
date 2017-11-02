import psycopg2,StringIO
from rdkit import Chem
from rdkit.Chem import Draw,AllChem


def get_conn():
    conn = psycopg2.connect(database='dsi',port=5433,host='52.91.71.182',user="postgres")
    return conn

def get_mol_list(results):
    mols = [Chem.MolFromSmiles(x[1]) for x in results]
    legends = [" ".join([str(x[0]),str(x[2])]) for x in results]
    return mols,legends

def ret_png(results):
    mols,legends = get_mol_list(results)
    img = Draw.MolsToGridImage(mols,legends=legends,molsPerRow=3)
    output = StringIO.StringIO()
    img.save(output, format="PNG")
    return output.getvalue()


def ret_svg(results):
    mols,legends = get_mol_list(results)
    img = Draw.MolsToGridImage(mols,legends=legends,molsPerRow=3,useSVG=True)
    return img


def get_graph_mol_list(results):
    mols = []
    legends = []
    for key in results:
        mol = Chem.MolFromSmiles(key.split("_")[0])
        mols.append(mol)
        legends.append("NODE")
        for m in results[key]:
            mols.append(Chem.MolFromSmiles(m))
            legends.append(key)
    return mols,legends

def ret_graph_svg(results):
    mols,legends = get_graph_mol_list(results)
    img = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=3, useSVG=True)
    return img

def draw_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "None Mol"
    AllChem.Compute2DCoords(mol)
    Chem.Kekulize(mol)
    drawer = Draw.MolDraw2DSVG(200,200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace('svg:','')