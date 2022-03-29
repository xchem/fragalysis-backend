from rest_framework import serializers

# Import standard models
from .models import Project, MculeQuote, Batch, Target, Method, Reaction, Product, Reactant, CatalogEntry, AnalyseAction

# Import action models
from .models import (
    AddAction,
    ExtractAction,
    FilterAction,
    QuenchAction,
    SetTemperatureAction,
    StirAction,
)

# Import OT session models
from .models import OTSession, Deck, Pipette, TipRack, Plate, Well, OTProtocol, OTBatchProtocol, CompoundOrder, OTScript

class CatalogEntrySerializer(serializers.ModelSerializer):
    class Meta:
        model = CatalogEntry
        fields = "__all__"

class ReactantSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reactant
        fields = "__all__"

class ReactantSerializerAll(serializers.ModelSerializer):
    catalogentries = CatalogEntrySerializer(many=True, read_only=True)

    class Meta:
        model = Reactant
        fields = "__all__"

class ProductSerializer(serializers.ModelSerializer):
    class Meta:
        model = Product
        fields = "__all__"

class ReactionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reaction
        fields = "__all__"

class ReactionSerializerAll(serializers.ModelSerializer):
    products = ProductSerializer(many=True, read_only=True)
    reactants = ReactantSerializerAll(many=True, read_only=True)

    class Meta:
        model = Reaction
        fields = "__all__"

class MethodSerializer(serializers.ModelSerializer):
    class Meta:
        model = Method
        fields = "__all__"

class MethodSerializerAll(serializers.ModelSerializer):
    reactions = ReactionSerializerAll(many=True, read_only=True)

    class Meta:
        model = Method
        fields = "__all__"

class TargetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Target
        fields = "__all__"

class TargetSerializerAll(serializers.ModelSerializer):
    methods = MethodSerializerAll(many=True, read_only=True)
    catalogentries = CatalogEntrySerializer(many=True, read_only=True)

    class Meta:
        model = Target
        fields = "__all__"

class BatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = Batch
        fields = "__all__"


class BatchSerializerAll(serializers.ModelSerializer):
    targets = TargetSerializerAll(many=True, read_only=True)

    class Meta:
        model = Batch
        fields = "__all__"


class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = "__all__"


class ProjectSerializerAll(serializers.ModelSerializer):
    batches = BatchSerializer(many=True, read_only=True)

    class Meta:
        model = Project
        fields = "__all__"



class MculeQuoteSerializer(serializers.ModelSerializer):
    class Meta:
        model = MculeQuote
        fields = "__all__"

# Action models here
class AnalyseActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalyseAction
        fields = "__all__"


class AddActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = AddAction
        fields = "__all__"


class ExtractActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExtractAction
        fields = "__all__"


class FilterActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = FilterAction
        fields = "__all__"


class QuenchActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = QuenchAction
        fields = "__all__"


class SetTemperatureActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = SetTemperatureAction
        fields = "__all__"


class StirActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = StirAction
        fields = "__all__"


# OT Session serializers
class DeckSerializer(serializers.ModelSerializer):
    class Meta:
        model = Deck
        fields = "__all__"


class PipetteSerializer(serializers.ModelSerializer):
    class Meta:
        model = Pipette
        fields = "__all__"


class TipRackSerializer(serializers.ModelSerializer):
    class Meta:
        model = TipRack
        fields = "__all__"


class PlateSerializer(serializers.ModelSerializer):
    class Meta:
        model = Plate
        fields = "__all__"


class WellSerializer(serializers.ModelSerializer):
    class Meta:
        model = Well
        fields = "__all__"


class CompoundOrderSerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundOrder
        fields = "__all__"


class OTScriptSerializer(serializers.ModelSerializer):
    class Meta:
        model = OTScript
        fields = "__all__"

class OTSessionSerializer(serializers.ModelSerializer):
    otscripts = OTScriptSerializer(many=True, read_only=True)
    compoundorders = CompoundOrderSerializer(many=True, read_only=True)

    class Meta:
        model = OTSession
        fields = "__all__"

class OTBatchProtocolSerializer(serializers.ModelSerializer):
    class Meta:
        model = OTBatchProtocol
        fields = "__all__"

class OTProtocolSerializer(serializers.ModelSerializer):
    class Meta:
        model = OTProtocol
        fields = "__all__"
