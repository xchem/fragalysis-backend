"""Targeted page-cache invalidation by view key_prefix.

Each cached view in viewer/views.py passes a unique key_prefix to cache_page.
clear_view_cache(*prefixes) drops the entries matching those prefixes — using
Redis SCAN/DELETE when the active cache backend is Redis, or clearing the
whole alias for backends without pattern delete.
"""
import logging

from django.conf import settings
from django.core.cache import caches
from django.core.cache.backends.redis import RedisCache

logger = logging.getLogger(__name__)


def clear_view_cache(*prefixes: str) -> None:
    """Invalidate cached pages whose key_prefix matches one of `prefixes`.

    With Redis, scans and deletes only the affected keys. With other backends
    (LocMemCache, DummyCache) there's no pattern delete, so the whole alias
    is cleared instead — coarser but correct.
    """
    if not prefixes:
        return
    cache = caches[settings.CACHE_MIDDLEWARE_ALIAS]
    if not isinstance(cache, RedisCache):
        cache.clear()
        return

    # pylint: disable=protected-access
    client = cache._cache.get_client(write=True)
    for prefix in prefixes:
        # cache_page stores per-user content under "...cache_page.{prefix}.{...}"
        # and the Vary lookup under "...cache_header.{prefix}.{...}"; clear both.
        for pattern in (f"*.cache_page.{prefix}.*", f"*.cache_header.{prefix}.*"):
            for key in client.scan_iter(match=pattern, count=500):
                client.delete(key)
