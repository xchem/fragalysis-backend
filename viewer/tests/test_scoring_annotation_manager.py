"""Regression test for the SiteObservationAnnotation manager (issue #988).

``SiteObservationAnnotationDataManager.get_queryset`` previously instantiated
``SiteObservationChoiceQueryset`` (a copy-paste), so
``SiteObservationAnnotation.filter_manager.filter_qs()`` - used directly as the
``SiteObservationAnnotationView`` queryset - built a query over the wrong model
(``SiteObservationChoice``). This pins the manager to its own model.
"""
from scoring.models import SiteObservationAnnotation, SiteObservationChoice


def test_annotation_filter_qs_uses_annotation_model():
    """filter_qs() must query SiteObservationAnnotation, not SiteObservationChoice."""
    assert (
        SiteObservationAnnotation.filter_manager.filter_qs().model
        is SiteObservationAnnotation
    )


def test_choice_filter_qs_still_uses_choice_model():
    """Guard the sibling manager so the fix does not swap the wrong way."""
    assert (
        SiteObservationChoice.filter_manager.filter_qs().model is SiteObservationChoice
    )
