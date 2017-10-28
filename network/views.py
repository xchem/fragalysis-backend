from django.http import HttpResponse
import json
from frag.network.query import get_picks

def index(request):
    if "smiles" in request.GET \
            and "num_picks" in request.GET:
        smiles = request.GET["smiles"]
        num_picks = request.GET["smiles"]
        out_dict = get_picks(smiles, num_picks)
        return HttpResponse(json.dumps(out_dict))
    else:
        return HttpResponse("Please insert SMILES")
