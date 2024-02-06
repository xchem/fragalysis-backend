from django.conf import settings
from django.http import HttpResponse
from frag.network.decorate import get_add_del_link
from frag.network.query import get_full_graph

from network.functions import order_structures


def full_graph(request):
    """Get the full graph for a molecule from an input smiles"""
    graph_choice = settings.NEO4J_QUERY
    graph_auth = settings.NEO4J_AUTH
    if "graph_choice" in request.GET:
        graph_choice = request.GET["graph_choice"]
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        out_dict = get_full_graph(smiles, graph_choice, graph_auth)
        decoration_list = get_add_del_link(smiles)
        if not out_dict:
            return HttpResponse("EMPTY RESULT SET")
        return HttpResponse(order_structures(out_dict, decoration_list))
    else:
        return HttpResponse("Please insert SMILES")
