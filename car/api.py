from pyexpat import model
from unittest.mock import DEFAULT
from django.db import DEFAULT_DB_ALIAS
from rest_framework import viewsets
from django.http import JsonResponse
from django.core.files.base import ContentFile
from django.contrib.admin.utils import NestedObjects
from django.db.models.fields.related import ForeignKey      

# Import standard models
from .models import Project, MculeQuote, Batch, Target, Method, Reaction, Product, AnalyseAction

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
    MculeQuoteSerializer,
    BatchSerializer,
    TargetSerializer,
    MethodSerializer,
    ReactionSerializer,
    ProductSerializer,
    AnalyseActionSerializer,
)

# Import action serializers
from .serializers import (
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

from collections import OrderedDict

def duplicatemethod(method: Method, new_target: Target):
    
    related_reaction_objs = method.reaction_set.all()

    # Duplicate method 
    method.pk = None
    method.target_id = new_target
    method.save()

    for reaction_obj in related_reaction_objs:
        reaction_obj.pk = None
        reaction_obj.method_id = method
        reaction_obj.save()

        # Reaction can only have one product
        product_obj = reaction_obj.product_set.all()[0]
        product_obj.pk = None
        product_obj.method_id = reaction_obj
        product_obj.save()

        # Reaction can have multiple add actions
        related_addaction_objs = reaction_obj.addaction_set.all()
        for addaction_obj in related_addaction_objs:
            addaction_obj.pk = None
            addaction_obj.method_id = reaction_obj
            addaction_obj.save()

        # Reaction can have multiple stir actions
        related_stiraction_objs = reaction_obj.stiraction_set.all()
        for stiraction_obj in related_stiraction_objs:
            stiraction_obj.pk = None
            stiraction_obj.method_id = reaction_obj
            stiraction_obj.save()

    return method.id
        
class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    serializer_class = ProjectSerializer


class MculeQuoteViewSet(viewsets.ModelViewSet):
    queryset = MculeQuote.objects.all()
    serializer_class = MculeQuoteSerializer

class BatchViewSet(viewsets.ModelViewSet):
    queryset = Batch.objects.all()
    serializer_class = BatchSerializer
    filterset_fields  = ["project_id"]


class TargetViewSet(viewsets.ModelViewSet):
    queryset = Target.objects.all()
    serializer_class = TargetSerializer
    filterset_fields  = ["batch_id"]


class MethodViewSet(viewsets.ModelViewSet):
    queryset = Method.objects.all()
    serializer_class = MethodSerializer
    filterset_fields = ["target_id", "nosteps"]


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    serializer_class = ReactionSerializer
    filterset_fields = {"method_id":["exact"] ,"successrate": ["gte", "lte"]}


class ProductViewSet(viewsets.ModelViewSet):
    queryset = Product.objects.all()
    serializer_class = ProductSerializer
    filterset_fields  = ["reaction_id"]


class AnalyseActionViewSet(viewsets.ModelViewSet):
    queryset = AnalyseAction.objects.all()
    serializer_class = AnalyseActionSerializer
    filterset_fields = ["reaction_id"]


# IBM viewsets here
class AddActionViewSet(viewsets.ModelViewSet):
    queryset = AddAction.objects.all()
    serializer_class = AddActionSerializer
    filterset_fields  = ["reaction_id"]

    def get_patch_object(self, pk):
        return AddAction.objects.get(pk=pk)

    # @action(methods=["PATCH"], detail=True)
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
    filterset_fields = ["project_id"]


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
    filterset_fields = ["project_id"]


class OTScriptViewSet(viewsets.ModelViewSet):
    queryset = OTScript.objects.all()
    serializer_class = OTScriptSerializer
    filterset_fields = ["otsession_id"]
