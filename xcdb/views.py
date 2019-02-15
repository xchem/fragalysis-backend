from xchem_db.models import (
    Crystal,
    DataProcessing,
    Dimple,
    Lab,
    Refinement,
    PanddaAnalysis,
    PanddaRun,
    PanddaSite,
    PanddaEvent,
    ProasisOut,
    PanddaEventStats,
)
from xchem_db.serializers import (
    CrystalSerializer,
    DataProcessingSerializer,
    DimpleSerializer,
    LabSerializer,
    RefinementSerializer,
    PanddaAnalysisSerializer,
    PanddaRunSerializer,
    PanddaSiteSerializer,
    PanddaEventSerializer,
    ProasisOutSerializer,
    PanddaEventStatsSerializer,
    FragspectCrystalView,
)

from rest_framework.views import APIView
from rest_framework.response import Response

from api.security import ISpyBSafeQuerySet


class CrystalView(ISpyBSafeQuerySet):
    queryset = Crystal.objects.filter()
    filter_permissions = "visit__proposal"
    serializer_class = CrystalSerializer
    filter_fields = (
        "crystal_name",
        "target__target_name",
        "compound__smiles",
        "visit__filename",
        "visit__proposal__proposal",
        "visit__visit",
    )


class DataProcessingView(ISpyBSafeQuerySet):
    queryset = DataProcessing.objects.filter()
    filter_permissions = "crystal_name__visit__proposal"
    serializer_class = DataProcessingSerializer
    filter_fields = (
        "crystal_name__crystal_name",
        "crystal_name__target__target_name",
        "crystal_name__compound__smiles",
        "crystal_name__visit__filename",
        "crystal_name__visit__proposal__proposal",
        "crystal_name__visit__visit",
    )


class DimpleView(ISpyBSafeQuerySet):
    queryset = Dimple.objects.filter()
    filter_permissions = "crystal_name__visit__proposal"
    serializer_class = DimpleSerializer
    filter_fields = (
        "crystal_name__crystal_name",
        "crystal_name__target__target_name",
        "crystal_name__compound__smiles",
        "crystal_name__visit__filename",
        "crystal_name__visit__proposal__proposal",
        "crystal_name__visit__visit",
        "reference__reference_pdb",
    )


class LabView(ISpyBSafeQuerySet):
    queryset = Lab.objects.filter()
    filter_permissions = "crystal_name__visit__proposal"
    serializer_class = LabSerializer
    filter_fields = (
        "crystal_name__crystal_name",
        "crystal_name__target__target_name",
        "crystal_name__compound__smiles",
        "crystal_name__visit__filename",
        "crystal_name__visit__proposal__proposal",
        "crystal_name__visit__visit",
        "data_collection_visit",
        "library_name",
        "library_plate",
    )


class RefinementView(ISpyBSafeQuerySet):
    queryset = Refinement.objects.filter()
    filter_permissions = "crystal_name__visit__proposal"
    serializer_class = RefinementSerializer
    filter_fields = (
        "crystal_name__crystal_name",
        "crystal_name__target__target_name",
        "crystal_name__compound__smiles",
        "crystal_name__visit__filename",
        "crystal_name__visit__proposal__proposal",
        "crystal_name__visit__visit",
        "outcome",
    )


class PanddaEventView(ISpyBSafeQuerySet):
    queryset = PanddaEvent.objects.filter()
    filter_permissions = "crystal__visit__proposal"
    serializer_class = PanddaEventSerializer
    filter_fields = (
        "crystal__crystal_name",
        "crystal__target__target_name",
        "crystal__compound__smiles",
        "crystal__visit__filename",
        "crystal__visit__proposal__proposal",
        "crystal__visit__visit",
        # "pandda_run__pandda_analysis__pandda_dir",
        # "pandda_run__pandda_log",
        # "pandda_run__sites_file",
        # "pandda_run__events_file",
        # "pandda_run__input_dir",
        # "site__site",
        # "event",
        # "lig_id",
        # "pandda_event_map_native",
        # "pandda_model_pdb",
        # "pandda_input_mtz",
        # "pandda_input_pdb",
    )


class PanddaEventStatsView(ISpyBSafeQuerySet):
    queryset = PanddaEventStats.objects.filter()
    filter_permissions = 'event__crystal__visit__proposal'
    serializer_class = PanddaEventStatsSerializer
    filter_fields = (
        "event__crystal__crystal_name",
        "event__crystal__target__target_name",
        "event__crystal__compound__smiles",
        "event__crystal__visit__filename",
        "event__crystal__visit__proposal__proposal",
        "event__crystal__visit__visit",
    )


class ProasisOutView(ISpyBSafeQuerySet):
    queryset = ProasisOut.objects.filter()
    filter_permissions = "crystal__visit__proposal"
    serializer_class = ProasisOutSerializer
    filter_fields = (
        "crystal__crystal_name",
        "crystal__target__target_name",
        "crystal__compound__smiles",
        "crystal__visit__filename",
        "crystal__visit__proposal__proposal",
        "crystal__visit__visit",
        "proasis__strucid",
        "proasis__crystal_name__crystal_name",
        "proasis__crystal_name__target__target_name",
        "proasis__crystal_name__compound__smiles",
        "proasis__crystal_name__visit__filename",
        "proasis__crystal_name__visit__proposal__proposal",
        "proasis__crystal_name__visit__visit",
        "proasis__refinement__crystal_name__crystal_name",
        "proasis__refinement__crystal_name__target__target_name",
        "proasis__refinement__crystal_name__compound__smiles",
        "proasis__refinement__crystal_name__visit__filename",
        "proasis__refinement__crystal_name__visit__proposal__proposal",
        "proasis__refinement__crystal_name__visit__visit",
        "proasis__refinement__outcome",
        "root",
        "start",
    )


class FragspectCrystalView(ISpyBSafeQuerySet):
    queryset = Refinement.objects.filter()
    serializer_class = FragspectCrystalView
    filter_fields = ('crystal_name__target__target_name',)
