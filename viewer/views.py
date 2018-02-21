from django.http import HttpResponse
from network.functions import draw_mol
from viewer.functions import generate_data_for_smis
from django.shortcuts import render
from .models import ViewScene,Molecule,Protein
from uuid import uuid4
import json
from frag.network.decorate import get_3d_vects_for_mol
from frag.network.query import get_full_graph
from fragalysis.utils import get_token
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger


def get_mols_from_scene(scene):
    comps = json.loads(scene)['components']
    mol_pks = []
    for comp in comps:
        path = comp['file_path']
        if "/viewer/mol_from_pk/" in path:
            mol_pks.append(int(path.split("/viewer/mol_from_pk/")[-1]))
    return mol_pks

def display(request):
    # Define the proteins, targets and molecules to be displayed / have options for displaying
    scene_id = 0
    if "target_title" in request.GET:
        target_title = request.GET["target_title"]
        mols = Molecule.objects.filter(prot_id__target_id__title=target_title)
    elif "scene_id" in request.GET:
        scene_id = int(request.GET["scene_id"])
        vs = ViewScene.objects.get(pk=scene_id)
        mol_pks = get_mols_from_scene(vs.scene)
        mols = Molecule.objects.filter(id__in=mol_pks)
    token = get_token(request)
    return render(request, 'viewer/display.html', {"token": token, "mols": mols, "scene_id": scene_id})

def inspect(request, target_pk):
    num_per_page = 5
    page = request.GET.get('page', 1)
    mol_pks = [x.pk for x in Molecule.objects.filter(prot_id__target_id__pk=target_pk) if x.prot_id.map_info.name != '']
    mols = Molecule.objects.filter(pk__in=mol_pks)
    paginator = Paginator(mols, num_per_page)
    try:
        mols = paginator.page(page)
    except PageNotAnInteger:
        mols = paginator.page(1)
    except EmptyPage:
        mols = paginator.page(paginator.num_pages)
    token = get_token(request)
    return render(request, 'viewer/inspect.html', {"token": token, "mols": mols})

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

def mol_from_pk(request, pk):
    sdf_info = Molecule.objects.get(pk=pk).sdf_info
    return  HttpResponse(sdf_info)


def map_from_pk(request,pk):
    prot = Protein.objects.get(pk=pk)
    return HttpResponse(prot.map_info)

def prot_from_pk(request,pk):
    pdb_info = open(Protein.objects.get(pk=pk).pdb_info.path).read()
    return  HttpResponse(pdb_info)


def post_view(request):
    """
    Post the view for a given scene
    :param request:
    :return:
    """
    new_view = ViewScene()
    new_view.title = request.POST["title"]
    new_view.scene = request.POST["scene"]
    new_view.uuid = str(uuid4())
    new_view.save()
    url = request.build_absolute_uri("?scene_id="+str(new_view.pk)).replace("/post_view/","/display/")
    return HttpResponse(url)

def get_view(request, pk):
    """
    Now Get the view for a given UUID
    :param request:
    :return:
    """
    this_view = ViewScene.objects.get(pk=pk)
    return HttpResponse(json.dumps({"title": this_view.title,
                         "scene": this_view.scene}))

def get_vects_from_pk(request, pk):
    sdf_info = str(Molecule.objects.get(pk=pk).sdf_info)
    out_data = get_3d_vects_for_mol(sdf_info)
    return HttpResponse(json.dumps(out_data))

def get_graph_from_pk(request,pk):
    smiles = str(Molecule.objects.get(pk=pk).smiles)
    out_data = get_full_graph(smiles)
    return HttpResponse(json.dumps(out_data))


def get_mols(request):
    mol_smiles = request.POST["mol_smiles"]
    if "target" in request.POST:
        target_id = request.POST["target"]
    out_data = generate_data_for_smis(mol_smiles)
    return HttpResponse(json.dumps(out_data))
