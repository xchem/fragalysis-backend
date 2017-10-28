from django.http import HttpResponse
from frag.network import

def index(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        return HttpResponse("Hello, world. You're at the polls index.")
    else:
        return HttpResponse("Please insert SMILES")
