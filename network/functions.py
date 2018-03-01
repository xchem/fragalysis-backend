import psycopg2
try:
    import StringIO
except:
    from io import StringIO
from rdkit import Chem
from rdkit.Chem import Draw,AllChem
import xml.etree.ElementTree as ET


def transparentsvg(svg):
    # Make the white background transparent
    tree = ET.fromstring(svg)
    rect = tree.find('rect')
    rect.set('style', rect.get('style').replace('#FFFFFF', 'none'))
    # Recover some missing attributes for correct browser rendering
    tree.set('version', '1.1')
    tree.set('xmlns', 'http://www.w3.org/2000/svg')
    tree.set('xmlns:rdkit', 'http://www.rdkit.org/xml')
    tree.set('xmlns:xlink', 'http://www.w3.org/1999/xlink')
    return '<?xml version="1.0" encoding="UTF-8"?>' + ET.tostring(tree).strip()

def get_conn():
    conn = psycopg2.connect(database='dsi',port=5432,host='cartridge',user="postgres")
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

def draw_mol(smiles,height=200,width=200):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return "None Mol"
    AllChem.Compute2DCoords(mol)
    Chem.Kekulize(mol)
    if not height:
        height=200
    if not width:
        width=200
    drawer = Draw.MolDraw2DSVG(height,width)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return transparentsvg(drawer.GetDrawingText().replace('svg:',''))

def add_keys(out_d,depth,type):
    if depth not in out_d:
        out_d[depth] = {}
    if type not in out_d[depth]:
        out_d[depth][type] = {}
    return out_d

def add_empty(out_d,tot_position,this_list,depth,this_type,annotation):
    for item in this_list:
        if item not in tot_position:
            out_d = add_keys(out_d, depth, this_type)
            out_d[depth][this_type][item] = {"smiles": [],"annotation":annotation}
    return out_d

def order_stuctures(results,decoration_list):
    """
    Order the data
    :param results:
    :return:
    """
    import json
    out_d = {}
    tot_position = []
    for key in results:
        depth = key.split("_")[1]
        this_type = key.split("_")[2]
        position = key.split("_")[0]
        tot_position.append(position)
        out_d = add_keys(out_d,depth,this_type)
        annotation = "BLANK"
        if position in decoration_list[0]:
            annotation = "ADD_DEC"
        if position in decoration_list[1]:
            annotation = "DEL_DEC"
        if position in decoration_list[2]:
            annotation = "LINK_DEC"
        out_d[depth][this_type][position] = {"smiles": results[key],"annotation":annotation}
    depth = -1
    out_d=add_empty(out_d,tot_position,decoration_list[0],depth,"ADDITION","ADD_MISS")
    out_d=add_empty(out_d, tot_position, decoration_list[1], depth, "DELETION","DEL_MISS")
    out_d=add_empty(out_d, tot_position, decoration_list[2], depth, "LINKER","LINK_MISS")
    return json.dumps(out_d)


def get_results(smiles,limit=100):
    conn = get_conn()
    curs = conn.cursor()
    curs.execute('select * from get_mfp2_neighbors(%s) limit ' + str(limit),
                 (smiles,))
    results = curs.fetchall()
    return results