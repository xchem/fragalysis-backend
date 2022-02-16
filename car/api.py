# Import standard models
from .models import Project, MculeQuote, Batch, Target, Method, Reaction, Product, AnalyseAction

# Import IBM models
from .models import (
    IBMAddAction,
    IBMCollectLayerAction,
    IBMConcentrateAction,
    IBMDegasAction,
    IBMDrySolidAction,
    IBMDrySolutionAction,
    IBMExtractAction,
    IBMFilterAction,
    IBMMakeSolutionAction,
    IBMPartitionAction,
    IBMpHAction,
    IBMPhaseSeparationAction,
    IBMQuenchAction,
    IBMRefluxAction,
    IBMSetTemperatureAction,
    IBMStirAction,
    IBMStoreAction,
    IBMWaitAction,
    IBMWashAction,
)

# Import OT Session models
from .models import OTSession, Deck, Pipette, TipRack, Plate, Well, CompoundOrder, OTScript

from rest_framework import viewsets
from django.http import JsonResponse
from django.core.files.base import ContentFile
from django.db.models import Q

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

# Import IBM serializers
from .serializers import (
    IBMAddActionSerializer,
    IBMCollectLayerActionSerializer,
    IBMConcentrateActionSerializer,
    IBMDegasActionSerializer,
    IBMDrySolidActionSerializer,
    IBMDrySolutionActionSerializer,
    IBMExtractActionSerializer,
    IBMFilterActionSerializer,
    IBMMakeSolutionActionSerializer,
    IBMPartitionActionSerializer,
    IBMpHActionSerializer,
    IBMPhaseSeparationActionSerializer,
    IBMQuenchActionSerializer,
    IBMRefluxActionSerializer,
    IBMSetTemperatureActionSerializer,
    IBMStirActionSerializer,
    IBMStoreActionSerializer,
    IBMWaitActionSerializer,
    IBMWashActionSerializer,
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


class GroupByStepsViewSet(viewsets.ModelViewSet):
    """ "
    Viewset to filter methods that have a given no of steps in a particular project
    """

    serializer_class = MethodSerializer

    def get_queryset(self):
        projectid = self.request.GET.get("projectid")
        nosteps = self.request.GET.get("nosteps")
        targetqueryset = Target.objects.filter(project_id=projectid).order_by("id")
        targetids = [target.id for target in targetqueryset]
        methodsqueryset = Method.objects.filter(Q(target_id__in=targetids) & Q(nosteps=nosteps))
        return methodsqueryset


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
class IBMAddActionViewSet(viewsets.ModelViewSet):
    queryset = IBMAddAction.objects.all()
    serializer_class = IBMAddActionSerializer
    filterset_fields  = ["reaction_id"]

    def get_patch_object(self, pk):
        return IBMAddAction.objects.get(pk=pk)

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
            serialized_data = IBMAddActionSerializer(addaction).data
            if serialized_data:
                return JsonResponse(data=serialized_data)
            else:
                return JsonResponse(data="wrong parameters")
        else:
            serializer = IBMAddActionSerializer(addaction, data=request.data, partial=True)
            if serializer.is_valid():
                serializer.save()
                return JsonResponse(data=serializer.data)
            else:
                return JsonResponse(data="wrong parameters")


class IBMCollectLayerActionViewSet(viewsets.ModelViewSet):
    queryset = IBMCollectLayerAction.objects.all()
    serializer_class = IBMCollectLayerActionSerializer
    filterset_fields = ["reaction_id"]


class IBMConcentrateActionViewSet(viewsets.ModelViewSet):
    queryset = IBMConcentrateAction.objects.all()
    serializer_class = IBMConcentrateActionSerializer
    filterset_fields = ["reaction_id"]


class IBMDegasActionViewSet(viewsets.ModelViewSet):
    queryset = IBMDegasAction.objects.all()
    serializer_class = IBMDegasActionSerializer
    filterset_fields = ["reaction_id"]


class IBMDrySolidActionViewSet(viewsets.ModelViewSet):
    queryset = IBMDrySolidAction.objects.all()
    serializer_class = IBMDrySolidActionSerializer
    filterset_fields = ["reaction_id"]


class IBMDrySolutionActionViewSet(viewsets.ModelViewSet):
    queryset = IBMDrySolutionAction.objects.all()
    serializer_class = IBMDrySolutionActionSerializer
    filterset_fields = ["reaction_id"]



class IBMExtractActionViewSet(viewsets.ModelViewSet):
    queryset = IBMExtractAction.objects.all()
    serializer_class = IBMExtractActionSerializer
    filterset_fields = ["reaction_id"]



class IBMFilterActionViewSet(viewsets.ModelViewSet):
    queryset = IBMFilterAction.objects.all()
    serializer_class = IBMFilterActionSerializer
    filterset_fields = ["reaction_id"]



class IBMMakeSolutionActionViewSet(viewsets.ModelViewSet):
    queryset = IBMMakeSolutionAction.objects.all()
    serializer_class = IBMMakeSolutionActionSerializer
    filterset_fields = ["reaction_id"]



class IBMPartitionActionViewSet(viewsets.ModelViewSet):
    queryset = IBMPartitionAction.objects.all()
    serializer_class = IBMPartitionActionSerializer
    filterset_fields = ["reaction_id"]



class IBMpHActionViewSet(viewsets.ModelViewSet):
    queryset = IBMpHAction.objects.all()
    serializer_class = IBMpHActionSerializer
    filterset_fields = ["reaction_id"]



class IBMPhaseSeparationActionViewSet(viewsets.ModelViewSet):
    queryset = IBMPhaseSeparationAction.objects.all()
    serializer_class = IBMPhaseSeparationActionSerializer
    filterset_fields = ["reaction_id"]


class IBMQuenchActionViewSet(viewsets.ModelViewSet):
    queryset = IBMQuenchAction.objects.all()
    serializer_class = IBMQuenchActionSerializer
    filterset_fields = ["reaction_id"]

class IBMRefluxActionViewSet(viewsets.ModelViewSet):
    queryset = IBMRefluxAction.objects.all()
    serializer_class = IBMRefluxActionSerializer
    filterset_fields = ["reaction_id"]


class IBMSetTemperatureActionViewSet(viewsets.ModelViewSet):
    queryset = IBMSetTemperatureAction.objects.all()
    serializer_class = IBMSetTemperatureActionSerializer
    filterset_fields = ["reaction_idd"]


class IBMStirActionViewSet(viewsets.ModelViewSet):
    queryset = IBMStirAction.objects.all()
    serializer_class = IBMStirActionSerializer
    filterset_fields = ["reaction_id"]


class IBMStoreActionViewSet(viewsets.ModelViewSet):
    queryset = IBMStoreAction.objects.all()
    serializer_class = IBMStoreActionSerializer
    filterset_fields = ["reaction_id"]


class IBMWaitActionViewSet(viewsets.ModelViewSet):
    queryset = IBMWaitAction.objects.all()
    serializer_class = IBMWaitActionSerializer
    filterset_fields = ["reaction_id"]


class IBMWashActionViewSet(viewsets.ModelViewSet):
    queryset = IBMWashAction.objects.all()
    serializer_class = IBMWashActionSerializer
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
