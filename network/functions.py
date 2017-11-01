import psycopg2,StringIO
from rdkit import Chem
from rdkit.Chem import Draw
def get_conn():
    conn = psycopg2.connect(database='dsi',port=5433,host='52.91.71.182',user="postgres")
    return conn



def get_mol_list(results):
    mols = [Chem.MolFromSmiles(x[1]) for x in results]
    legends = [" ".join([x[0],x[2]]) for x in results]
    return mols,legends

def ret_png(results):
    mols,legends = get_mol_list(results)
    img = Draw.MolsToGridImage(mols,legends=legends,molsPerRow=6)
    output = StringIO.StringIO()
    img.save(output, format="PNG")
    return output.getvalue()


def ret_svg(results):
    mols,legends = get_mol_list(results)
    img = Draw.MolsToGridImage(mols,legends=legends,molsPerRow=6,useSVG=True)
    output = StringIO.StringIO()
    img.save(output, format="SVG")
    return output.getvalue()