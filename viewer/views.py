from django.http import HttpResponse
import json
from network.functions import draw_mol
from django.shortcuts import render

def display(request):
    return render(request, 'viewer/display.html', {})

def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        return  HttpResponse(draw_mol(smiles))
    else:
        return HttpResponse("Please insert SMILES")