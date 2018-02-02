from django.http import HttpResponse
from network.functions import draw_mol
from django.shortcuts import render
from django.http import JsonResponse
from django.views import View
from .models import Protein,ViewScene
from .forms import PDBForm
from uuid import uuid4

def display(request):
    # Define the proteins, targets and molecules to be displayed / have options for displaying
    return render(request, 'viewer/display.html', {})

def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        return  HttpResponse(draw_mol(smiles))
    else:
        return HttpResponse("Please insert SMILES")


def post_view(request):
    """
    Post the view for a given scene
    :param request:
    :return:
    """
    new_view = ViewScene()
    new_view.title = request.POST["data"]["title"]
    new_view.scene = request.POST["data"]["title"]
    new_view.uuid = str(uuid4())
    new_view.save()
    return HttpResponse(new_view.uuid)

def get_view(request):
    """
    Now Get the view for a given UUID
    :param request:
    :return:
    """
    uuid = request.GET["uuid"]
    this_view = ViewScene.objects.get(uuid=uuid)
    return HttpResponse({"title": this_view.title,
                         "scene": this_view.scene})


# POST AND GETS FOR ALL THE OTHER FILE TYPES.

