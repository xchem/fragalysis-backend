from xchem_db.models import (
    Crystal,
    DataProcessing,
    Dimple,
    Lab,
    PanddaEvent,
    PanddaEventStats,
    ProasisOut,
    Refinement,
)
from xchem_db.serializers import (
    CrystalSerializer,
    DataProcessingSerializer,
    DimpleSerializer,
    FragspectCrystalSerializer,
    LabSerializer,
    PanddaEventSerializer,
    PanddaEventStatsSerializer,
    ProasisOutSerializer,
    RefinementSerializer,
)

from api.security import ISPyBSafeQuerySet


class CrystalView(ISPyBSafeQuerySet):
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


class DataProcessingView(ISPyBSafeQuerySet):
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


class DimpleView(ISPyBSafeQuerySet):
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


class LabView(ISPyBSafeQuerySet):
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


class RefinementView(ISPyBSafeQuerySet):
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


class PanddaEventView(ISPyBSafeQuerySet):
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


class PanddaEventStatsView(ISPyBSafeQuerySet):
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


class ProasisOutView(ISPyBSafeQuerySet):
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


class FragspectCrystalView(ISPyBSafeQuerySet):
    queryset = PanddaEvent.objects.filter().prefetch_related(
        "crystal__target",
        "crystal__compound",
        "crystal",
        "site",
        "refinement",
        "data_proc",
    )
    serializer_class = FragspectCrystalSerializer
    filter_fields = {"crystal__target__target_name": ["iexact"]}
    filter_permissions = "crystal__visit__proposal"
