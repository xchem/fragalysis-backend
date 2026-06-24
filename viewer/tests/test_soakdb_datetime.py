"""Unit tests for ``TargetLoader._soakdb_datetime`` timezone handling.

SoakDB timestamps carry no timezone, so ``dateutil.parser.parse`` returns a
*naive* datetime. With Django's ``USE_TZ`` active, storing a naive datetime
emits a ``RuntimeWarning`` on every ``Experiment`` datetime field - tens of
thousands of warnings per target load (visible in the integration test's Celery
worker log). The loader must therefore return timezone-aware datetimes so those
warnings never fire.

The datetime path touches no instance state, so each test builds the loader with
``object.__new__`` to skip the filesystem-heavy ``__init__``.
"""

# The method under test is intentionally "protected"; testing it directly is the
# point here, so silence pylint's protected-access for the module.
# pylint: disable=protected-access
from datetime import datetime

from django.utils import timezone

from viewer.target_loader import TargetLoader


def _loader() -> TargetLoader:
    """A TargetLoader whose ``__init__`` is bypassed (datetime path needs none)."""
    return object.__new__(TargetLoader)


def test_soakdb_datetime_makes_naive_value_aware():
    """A bare SoakDB timestamp (no offset) is returned timezone-aware."""
    result = _loader()._soakdb_datetime(
        {"DataCollectionDate": "2025-09-20 12:25:41"}, "DataCollectionDate"
    )
    assert isinstance(result, datetime)
    assert timezone.is_aware(result), "naive SoakDB datetime must be made aware"


def test_soakdb_datetime_empty_values_return_none():
    """Empty / literal-"None" fields still map to ``None`` (no datetime)."""
    loader = _loader()
    assert loader._soakdb_datetime({"RefinementDate": ""}, "RefinementDate") is None
    assert loader._soakdb_datetime({"RefinementDate": "None"}, "RefinementDate") is None
