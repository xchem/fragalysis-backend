from rest_framework import viewsets
from rest_framework.decorators import action
from django.http import JsonResponse
from django.core.files.base import ContentFile
from celery.result import AsyncResult
from viewer.tasks import check_services

from car.tasks import createOTScript  

# Import standard models
from .models import Project, MculeQuote, Batch, Target, Method, Reaction, Reactant, CatalogEntry, Product, AnalyseAction

# Import action models
from .models import (
    AddAction,
    ExtractAction,
    FilterAction,
    QuenchAction,
    SetTemperatureAction,
    StirAction,

)

# Import OT Session models
from .models import OTSession, Deck, Pipette, TipRack, Plate, Well, CompoundOrder, OTScript

# Import standard serializers
from .serializers import (
    ProjectSerializer,
    ProjectSerializerAll,
    MculeQuoteSerializer,
    BatchSerializer,
    BatchSerializerAll,
    TargetSerializer,
    TargetSerializerAll,
    MethodSerializer,
    MethodSerializerAll,
    ReactionSerializer,
    ReactionSerializerAll,
    ProductSerializer,
    ReactantSerializer,
    ReactantSerializerAll,
    CatalogEntrySerializer,
)

# Import action serializers
from .serializers import (
    AnalyseActionSerializer,
    AddActionSerializer,
    ExtractActionSerializer,
    FilterActionSerializer,  
    QuenchActionSerializer,
    SetTemperatureActionSerializer,
    StirActionSerializer,
    
)

# Import OT Session serializers
from .serializers import (
    OTSessionSerializer,
    DeckSerializer,
    PipetteSerializer,
    TipRackSerializer,
    PlateSerializer,
    WellSerializer,
    CompoundOrderSerializer,
    OTScriptSerializer,
)

from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage

from .utils import createSVGString

def duplicatetarget(target_obj: Target, fk_obj: Batch):
    related_catalogentry_queryset = target_obj.catalogentries.all()
    
    target_obj.image = ContentFile(target_obj.image.read(), name=target_obj.image.name)
    target_obj.pk = None
    target_obj.batch_id = fk_obj
    target_obj.save()

    for catalogentry_obj in related_catalogentry_queryset:
        catalogentry_obj.pk = None
        catalogentry_obj.target_id = target_obj
        catalogentry_obj.save()

    return target_obj

def duplicatemethod(method_obj: Method, fk_obj: Target):
    related_reaction_queryset = method_obj.reactions.all()

    # Duplicate method before cloning the related children reaction objs
    method_obj.pk = None
    method_obj.target_id = fk_obj
    method_obj.save()

    for reaction_obj in related_reaction_queryset:
        product_obj = reaction_obj.products.all()[0]
        related_addaction_objs = reaction_obj.addactions.all()
        related_stiraction_objs = reaction_obj.stiractions.all()
        related_reactant_objs = reaction_obj.reactants.all()
  
        reaction_obj.reactionimage = ContentFile(reaction_obj.reactionimage.read(), 
                                                name=reaction_obj.reactionimage.name)
        reaction_obj.pk = None
        reaction_obj.method_id = method_obj
        reaction_obj.save()

        product_obj.image = ContentFile(product_obj.image.read(), 
                                        name=product_obj.image.name)
        product_obj.pk = None
        product_obj.reaction_id = reaction_obj
        product_obj.save()

        for addaction_obj in related_addaction_objs:
            addaction_obj.pk = None
            addaction_obj.reaction_id = reaction_obj
            addaction_obj.save()

        for stiraction_obj in related_stiraction_objs:
            stiraction_obj.pk = None
            stiraction_obj.reaction_id = reaction_obj
            stiraction_obj.save()

        for reactant_obj in related_reactant_objs:
            related_catalogentry_objs = reactant_obj.catalogentries.all()
            reactant_obj.pk = None
            reactant_obj.reaction_id = reaction_obj
            reactant_obj.save()
            for catalog_obj in related_catalogentry_objs:
                catalog_obj.pk = None
                catalog_obj.reactant_id = reactant_obj
                catalog_obj.save()
        
        
class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    
    def get_serializer_class(self):
        fetchall = self.request.GET.get('fetchall', None)
        return ProjectSerializerAll if fetchall == "yes" else ProjectSerializer


class MculeQuoteViewSet(viewsets.ModelViewSet):
    queryset = MculeQuote.objects.all()
    serializer_class = MculeQuoteSerializer
    

class BatchViewSet(viewsets.ModelViewSet):
    queryset = Batch.objects.all()
    filterset_fields  = ["project_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get('fetchall', None)
        return BatchSerializerAll if fetchall == "yes" else BatchSerializer

    def createBatch(self, project_obj, batch_node_obj, batch_tag):
        batch_obj = Batch()
        batch_obj.project_id = project_obj
        batch_obj.batch_id = batch_node_obj
        batch_obj.batch_tag = batch_tag
        batch_obj.save()
        return batch_obj

    def create(self, request, **kwargs):
        method_ids = request.data["methodids"]
        batch_tag = request.data["batchtag"]
        try: 
            target_query_set = Target.objects.filter(methods__id__in=method_ids).distinct()
            batch_obj = target_query_set[0].batch_id
            project_obj = batch_obj.project_id
            batch_obj_new = self.createBatch(project_obj=project_obj, batch_node_obj=batch_obj, batch_tag=batch_tag)
            for target_obj in target_query_set:
                method_query_set_to_clone = Method.objects.filter(target_id=target_obj).filter(pk__in=method_ids)
                target_obj_clone = duplicatetarget(target_obj=target_obj, fk_obj=batch_obj_new)
                for method_obj in method_query_set_to_clone:
                    duplicatemethod(method_obj=method_obj, fk_obj=target_obj_clone)
            serialized_data = BatchSerializer(batch_obj_new).data
            if serialized_data:
                return JsonResponse(data=serialized_data)
            else:
                return JsonResponse(data="Something went wrong")
        except:
            return JsonResponse(data="Something went wrong")


class TargetViewSet(viewsets.ModelViewSet):
    queryset = Target.objects.all()
    filterset_fields  = ["batch_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get('fetchall', None)
        return TargetSerializerAll if fetchall == "yes" else TargetSerializer
    


class MethodViewSet(viewsets.ModelViewSet):
    queryset = Method.objects.all()
    filterset_fields = ["target_id", "nosteps"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get('fetchall', None)
        return MethodSerializerAll if fetchall == "yes" else MethodSerializer


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    filterset_fields = {"method_id":["exact"] ,"successrate": ["gte", "lte"]}

    def get_serializer_class(self):
        fetchall = self.request.GET.get('fetchall', None)
        return ReactionSerializerAll if fetchall == "yes" else ReactionSerializer

class ProductViewSet(viewsets.ModelViewSet):
    queryset = Product.objects.all()
    serializer_class = ProductSerializer
    filterset_fields  = ["reaction_id"]

class ReactantViewSet(viewsets.ModelViewSet):
    queryset = Reactant.objects.all()
    filterset_fields = ["reaction_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get('fetchall', None)
        return ReactantSerializerAll if fetchall == "yes" else ReactantSerializer

class CatalogEntryViewSet(viewsets.ModelViewSet):
    queryset = CatalogEntry.objects.all()
    serializer_class = CatalogEntrySerializer

# Action viewsets
class AnalyseActionViewSet(viewsets.ModelViewSet):
    queryset = AnalyseAction.objects.all()
    serializer_class = AnalyseActionSerializer
    filterset_fields = ["reaction_id"]


class AddActionViewSet(viewsets.ModelViewSet):
    queryset = AddAction.objects.all()
    serializer_class = AddActionSerializer
    filterset_fields  = ["reaction_id"]

    def get_patch_object(self, pk):
        return AddAction.objects.get(pk=pk)

    def partial_update(self, request, pk):
        addaction = self.get_patch_object(pk)

        if "materialsmiles" in request.data:
            materialsmiles = request.data["materialsmiles"]
            mol = Chem.MolFromSmiles(materialsmiles)
            molecular_weight = Descriptors.ExactMolWt(mol)
            add_svg_string = createSVGString(materialsmiles)
            # Delete previous image
            svg_fp = addaction.materialimage.path
            default_storage.delete(svg_fp)
            # Add new image
            add_svg_fn = default_storage.save(
                "addactionimages/{}.svg".format(materialsmiles),
                ContentFile(add_svg_string),
            )
            addaction.materialsmiles = materialsmiles
            addaction.molecularweight = molecular_weight
            addaction.materialimage = add_svg_fn
            addaction.save()
            serialized_data = AddActionSerializer(addaction).data
            if serialized_data:
                return JsonResponse(data=serialized_data)
            else:
                return JsonResponse(data="wrong parameters")
        else:
            serializer = AddActionSerializer(addaction, data=request.data, partial=True)
            if serializer.is_valid():
                serializer.save()
                return JsonResponse(data=serializer.data)
            else:
                return JsonResponse(data="wrong parameters")



class ExtractActionViewSet(viewsets.ModelViewSet):
    queryset = ExtractAction.objects.all()
    serializer_class = ExtractActionSerializer
    filterset_fields = ["reaction_id"]


class FilterActionViewSet(viewsets.ModelViewSet):
    queryset = FilterAction.objects.all()
    serializer_class = FilterActionSerializer
    filterset_fields = ["reaction_id"]


class QuenchActionViewSet(viewsets.ModelViewSet):
    queryset = QuenchAction.objects.all()
    serializer_class = QuenchActionSerializer
    filterset_fields = ["reaction_id"]


class SetTemperatureActionViewSet(viewsets.ModelViewSet):
    queryset = SetTemperatureAction.objects.all()
    serializer_class = SetTemperatureActionSerializer
    filterset_fields = ["reaction_idd"]


class StirActionViewSet(viewsets.ModelViewSet):
    queryset = StirAction.objects.all()
    serializer_class = StirActionSerializer
    filterset_fields = ["reaction_id"]


# OT Session viewsets
class OTSessionViewSet(viewsets.ModelViewSet):
    queryset = OTSession.objects.all()
    serializer_class = OTSessionSerializer
    filterset_fields = ["batch_id"]

    @action(methods=['post'], detail=False)
    def createotprotocol(self, request, pk=None):
        check_services()
        # batch_ids = request.data["batchids"]
        batch_ids = request.POST.getlist("batchids")
        batch_ids = [int(batchid) for batchid in batch_ids]
        print(batch_ids)
        task = createOTScript.delay(batchids=batch_ids)
        data = {"task_id": task.id}
        return JsonResponse(data=data)
    
    @action(detail=False, methods=['get'])
    def gettaskstatus(self, request, pk=None):
        task_id = self.request.GET.get('task_id', None)
        if task_id:
            task = AsyncResult(task_id)            
            if task.status == "FAILURE":
                data = {"task_status": task.status, "traceback": str(task.traceback)}
                return JsonResponse(data)

            if task.status == "SUCCESS":
                result = task.get()
                data = {"task_status": task.status, "protocol_summary": result}
                return JsonResponse(data)
                
            if task.status == "PENDING":
                data = {"task_status": task.status}
                return JsonResponse(data)

class DeckViewSet(viewsets.ModelViewSet):
    queryset = Deck.objects.all()
    serializer_class = DeckSerializer
    filterset_fields = ["otsession_id"]


class PipetteViewSet(viewsets.ModelViewSet):
    queryset = Pipette.objects.all()
    serializer_class = PipetteSerializer
    filterset_fields = ["otsession_id"]


class TipRackViewSet(viewsets.ModelViewSet):
    queryset = TipRack.objects.all()
    serializer_class = TipRackSerializer
    filterset_fields = ["otsession_id"]


class PlateViewSet(viewsets.ModelViewSet):
    queryset = Plate.objects.all()
    serializer_class = PlateSerializer
    filterset_fields = ["otsession_id"]


class WellViewSet(viewsets.ModelViewSet):
    queryset = Well.objects.all()
    serializer_class = WellSerializer
    filterset_fields = ["otsession_id"]


class CompoundOrderViewSet(viewsets.ModelViewSet):
    queryset = CompoundOrder.objects.all()
    serializer_class = CompoundOrderSerializer
    filterset_fields = ["otsession_id"]


class OTScriptViewSet(viewsets.ModelViewSet):
    queryset = OTScript.objects.all()
    serializer_class = OTScriptSerializer
    filterset_fields = ["otsession_id"]