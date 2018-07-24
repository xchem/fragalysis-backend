from django.http import HttpResponse
from frag.network.query import get_full_graph, custom_query
from frag.network.decorate import get_add_del_link
from network.functions import order_stuctures


def full_graph(request):
    """
    Get the full graph for a molecule from an input smiles
    :param request:
    :return:
    """
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        out_dict = get_full_graph(smiles)
        decoration_list = get_add_del_link(smiles)
        if not out_dict:
            return HttpResponse("EMPTY RESULT SET")
        return HttpResponse(order_stuctures(out_dict, decoration_list))
    else:
        return HttpResponse("Please insert SMILES")


def custom_query(request):
    query = request.GET.get("query", None)
    if query:
        return HttpResponse(custom_query(query))
    else:
        return HttpResponse("Please insert Query")
