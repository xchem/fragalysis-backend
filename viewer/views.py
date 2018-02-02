from django.http import HttpResponse
from network.functions import draw_mol
from django.shortcuts import render
from .models import ViewScene,Molecule,Protein
from uuid import uuid4

def display(request):
    # Define the proteins, targets and molecules to be displayed / have options for displaying
    target_title = request.GET["target_title"]
    mols = Molecule.objects.filter(prot_id__target_id__title=target_title)
    return render(request, 'viewer/display.html', {"mols":mols})

def mol_view(request):
    if "smiles" in request.GET:
        smiles = request.GET["smiles"].rstrip(".svg")
        return  HttpResponse(draw_mol(smiles))
    else:
        return HttpResponse("Please insert SMILES")

def img_from_pk(request):
    if "pk" in request.GET:
        smiles = Molecule.objects.get(pk=request.GET["pk"]).smiles
        return  HttpResponse(draw_mol(smiles))
    else:
        return HttpResponse("Please insert PK")

def mol_from_pk(request):
    if "pk" in request.GET:
        sdf_info = Molecule.objects.get(pk=request.GET["pk"]).sdf_info
        return  HttpResponse(sdf_info)
    else:
        return HttpResponse("Please insert PK")

def prot_from_pk(request):
    if "pk" in request.GET:
        pdb_info = open(Protein.objects.get(pk=request.GET["pk"]).pdb_info.path).read()
        return  HttpResponse(pdb_info)
    else:
        return HttpResponse("Please insert PK")


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

