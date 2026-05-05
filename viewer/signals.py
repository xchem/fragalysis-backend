from django.conf import settings
from django.core.cache import caches
from django.db.models.signals import post_delete, post_save

from viewer.models import ComputedMolecule, Pose, SiteObservation, SiteObservationTag

_CACHED_MODELS = (
    SiteObservation,
    ComputedMolecule,
    Pose,
    SiteObservationTag,
)


def _clear_view_cache(**_kwargs):
    # Note: bulk_create/bulk_update/queryset.update/queryset.delete bypass these
    # signals and must clear the cache explicitly at their call site.
    caches[settings.CACHE_MIDDLEWARE_ALIAS].clear()


for _model in _CACHED_MODELS:
    post_save.connect(_clear_view_cache, sender=_model, weak=False)
    post_delete.connect(_clear_view_cache, sender=_model, weak=False)
