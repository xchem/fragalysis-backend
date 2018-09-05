from django.http import HttpResponse
from frag.network.query import get_full_graph
from frag.network.decorate import get_add_del_link
from network.functions import order_stuctures
import os


def full_graph(request):
    """
    Get the full graph for a molecule from an input smiles
    :param request:
    :return:
    """
    graph_choice = os.environ.get("NEO4J_QUERY", "neo4j")
    if "graph_choice" in request.GET:
        graph_choice = request.GET["graph_choice"]
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        out_dict = get_full_graph(smiles, graph_choice)
        decoration_list = get_add_del_link(smiles)
        if not out_dict:
            return HttpResponse("EMPTY RESULT SET")
        return HttpResponse(order_stuctures(out_dict, decoration_list))
    else:
        return HttpResponse("Please insert SMILES")
