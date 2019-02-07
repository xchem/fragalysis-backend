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
)

from api.security import ISpyBSafeQuerySet


class CrystalView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = Crystal.objects.filter()
    filter_permissions = "visit__number"
    serializer_class = CrystalSerializer
    filter_fields = (
        "crystal_name",
        "target__target_name",
        "compound__smiles",
        "visit__filename",
        "visit__proposal__proposal",
        "visit__visit",
    )


class DataProcessingView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = DataProcessing.objects.filter()
    filter_permissions = "crystal__visit__number"
    serializer_class = DataProcessingSerializer
    filter_fields = (
        "crystal_name__crystal_name",
        "crystal_name__target__target_name",
        "crystal_name__compound__smiles",
        "crystal_name__visit__filename",
        "crystal_name__visit__proposal__proposal",
        "crystal_name__visit__visit",
    )


class DimpleView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = Dimple.objects.filter()
    filter_permissions = "crystal__visit__number"
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


class LabView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = Lab.objects.filter()
    filter_permissions = "crystal__visit__number"
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


class RefinementView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = Refinement.objects.filter()
    filter_permissions = "crystal__visit__number"
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


class PanddaAnalysisView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = PanddaAnalysis.objects.filter()
    filter_permissions = "crystal__visit__number"
    serializer_class = PanddaAnalysisSerializer
    filter_fields = ("pandda_dir",)


class PanddaRunView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = PanddaRun.objects.filter()
    filter_permissions = "crystal__visit__number"
    serializer_class = PanddaRunSerializer
    filter_fields = (
        "input_dir",
        "pandda_analysis__pandda_dir",
        "pandda_log",
        "pandda_version",
        "sites_file",
        "events_file",
    )


class PanddaSiteView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = PanddaSite.objects.filter()
    filter_permissions = "crystal__visit__number"
    serializer_class = PanddaSiteSerializer
    filter_fields = (
        "pandda_run__pandda_analysis__pandda_dir",
        "pandda_run__pandda_log",
        "pandda_run__pandda_sites_file",
        "pandda_run__pandda_events_file",
        "pandda_run__input_dir",
        "site",
    )


class PanddaEventView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = PanddaEvent.objects.filter()
    filter_permissions = "crystal__visit__number"
    serializer_class = PanddaEventSerializer
    filter_fields = (
        "crystal__crystal_name",
        "crystal__target__target_name",
        "crystal__compound__smiles",
        "crystal__visit__filename",
        "crystal__visit__proposal__proposal",
        "crystal__visit__visit",
        "pandda_run__pandda_analysis__pandda_dir",
        "pandda_run__pandda_log",
        "pandda_run__pandda_sites_file",
        "pandda_run__pandda_events_file",
        "pandda_run__input_dir",
        "site__site",
        "event",
        "lig_id",
        "pandda_event_map_native",
        "pandda_model_pdb",
        "pandda_input_mtz",
        "pandda_input_pdb",
    )


class ProasisOutView(ISpyBSafeQuerySet(xcdb=True)):
    queryset = ProasisOut.objects.filter()
    filter_permissions = "crystal__visit__number"
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
