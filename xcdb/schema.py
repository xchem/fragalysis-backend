import graphene
from graphene_django.rest_framework.mutation import SerializerMutation

from xchem_db.serializers import (TargetSerializer, CompoundsSerializer, ReferenceSerializer, SoakdbFilesSerializer,
                                  CrystalSerializer, DataProcessingSerializer, DimpleSerializer, LabSerializer,
                                  RefinementSerializer, PanddaAnalysisSerializer,
                                  PanddaRunSerializer, PanddaSiteSerializer, PanddaEventSerializer,
                                  ProasisOutSerializer, FragspectCrystalSerializer)

relay = graphene.relay


class Target(SerializerMutation):
    serializer_class = TargetSerializer
    interfaces = (relay.Node,)


class Compounds(SerializerMutation):
    serializer_class = CompoundsSerializer
    interfaces = (relay.Node,)


class Reference(SerializerMutation):
    serializer_class = ReferenceSerializer
    interfaces = (relay.Node,)


class SoakdbFiles(SerializerMutation):
    serializer_class = SoakdbFilesSerializer
    interfaces = (relay.Node,)


class Crystal(SerializerMutation):
    serializer_class = CrystalSerializer
    interfaces = (relay.Node,)


class DataProcessing(SerializerMutation):
    serializer_class = DataProcessingSerializer
    interfaces = (relay.Node,)


class Dimple(SerializerMutation):
    serializer_class = DimpleSerializer
    interfaces = (relay.Node,)


class Lab(SerializerMutation):
    serializer_class = LabSerializer
    interfaces = (relay.Node,)


class Refinement(SerializerMutation):
    serializer_class = RefinementSerializer
    interfaces = (relay.Node,)


class PanddaAnalysis(SerializerMutation):
    serializer_class = PanddaAnalysisSerializer
    interfaces = (relay.Node,)


class PanddaRun(SerializerMutation):
    serializer_class = PanddaRunSerializer
    interfaces = (relay.Node,)


class PanddaSite(SerializerMutation):
    serializer_class = PanddaSiteSerializer
    interfaces = (relay.Node,)


class PanddaEvent(SerializerMutation):
    serializer_class = PanddaEventSerializer
    interfaces = (relay.Node,)


class ProasisOut(SerializerMutation):
    serializer_class = ProasisOutSerializer
    interfaces = (relay.Node,)

class Fragspect(SerializerMutation):
    serializer_class = FragspectCrystalView
    interfaces = (relay.Node,)


class Query(graphene.ObjectType):
    target = graphene.list(Target)
    compounds = graphene.list(Compounds)
    reference = graphene.list(Reference)
    soakdb_files = graphene.list(SoakdbFiles)
    crystal = graphene.list(Crystal)
    data_processing = graphene.list(DataProcessing)
    dimple = graphene.list(Dimple)
    lab = graphene.list(Lab)
    refinement = graphene.list(Refinement)
    pandda_analysis = graphene.list(PanddaAnalysis)
    pandda_run = graphene.list(PanddaRun)
    pandda_site = graphene.list(PanddaSite)
    pandda_event = graphene.list(PanddaEvent)
    proasis_out = graphene.list(ProasisOut)
    fragspect = graphene.list(Fragspect)


schema = graphene.Schema(query=Query)
