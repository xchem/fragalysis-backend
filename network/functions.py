import json


def add_keys(out_d, depth, this_type):
    if depth not in out_d:
        out_d[depth] = {}
    if this_type not in out_d[depth]:
        out_d[depth][this_type] = {}
    return out_d


def add_empty(out_d, tot_position, this_list, depth, this_type, annotation):
    for item in this_list:
        if item not in tot_position:
            out_d = add_keys(out_d, depth, this_type)
            out_d[depth][this_type][item] = {"smiles": [], "annotation": annotation}
    return out_d


def order_structures(results, decoration_list):
    """
    Order the data
    :param results:
    :return:
    """
    out_d = {}
    tot_position = []
    for key in results:
        depth = key.split("_")[1]
        this_type = key.split("_")[2]
        position = key.split("_")[0]
        tot_position.append(position)
        out_d = add_keys(out_d, depth, this_type)
        annotation = "BLANK"
        if position in decoration_list[0]:
            annotation = "ADD_DEC"
        if position in decoration_list[1]:
            annotation = "DEL_DEC"
        if position in decoration_list[2]:
            annotation = "LINK_DEC"
        out_d[depth][this_type][position] = {
            "smiles": results[key],
            "annotation": annotation,
        }
    depth = -1
    out_d = add_empty(
        out_d, tot_position, decoration_list[0], depth, "ADDITION", "ADD_MISS"
    )
    out_d = add_empty(
        out_d, tot_position, decoration_list[1], depth, "DELETION", "DEL_MISS"
    )
    out_d = add_empty(
        out_d, tot_position, decoration_list[2], depth, "LINKER", "LINK_MISS"
    )
    return json.dumps(out_d)
