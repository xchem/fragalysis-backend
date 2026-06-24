"""Tests for the page-cache invalidation helpers in ``viewer.cache`` (#991).

``clear_view_cache`` deletes only the keys matching the given view key_prefixes
when the backend is Redis (SCAN/DELETE), and falls back to clearing the whole
alias for backends without pattern delete. The Django cache backend is faked so
the key-pattern logic is asserted without a real Redis.
"""
from unittest import mock

from django.core.cache.backends.redis import RedisCache

from viewer import cache as cache_module


class _FakeCaches:
    """Stands in for ``django.core.cache.caches`` - always hands out one cache."""

    def __init__(self, backend):
        self._backend = backend

    def __getitem__(self, _alias):
        return self._backend


def _patch_caches(monkeypatch, backend):
    monkeypatch.setattr(cache_module, "caches", _FakeCaches(backend))


def test_clear_view_cache_no_prefixes_is_noop(monkeypatch):
    backend = mock.Mock()
    _patch_caches(monkeypatch, backend)

    cache_module.clear_view_cache()

    backend.clear.assert_not_called()


def test_clear_view_cache_non_redis_clears_whole_alias(monkeypatch):
    backend = mock.Mock()  # not a RedisCache instance
    _patch_caches(monkeypatch, backend)

    cache_module.clear_view_cache("alpha")

    backend.clear.assert_called_once_with()


def test_clear_view_cache_redis_scans_and_deletes(monkeypatch):
    client = mock.Mock()
    # Two keys returned per scan pattern.
    client.scan_iter.return_value = ["k1", "k2"]

    backend = mock.Mock(spec=RedisCache)
    backend._cache.get_client.return_value = client  # pylint: disable=protected-access
    _patch_caches(monkeypatch, backend)

    cache_module.clear_view_cache("alpha")

    # cache_page + cache_header patterns are both scanned for the prefix.
    patterns = [c.kwargs["match"] for c in client.scan_iter.call_args_list]
    assert "*.cache_page.alpha.*" in patterns
    assert "*.cache_header.alpha.*" in patterns
    # Every scanned key is deleted (2 keys x 2 patterns).
    assert client.delete.call_count == 4
    backend.clear.assert_not_called()


def test_clear_all_view_caches_clears_alias(monkeypatch):
    backend = mock.Mock()
    _patch_caches(monkeypatch, backend)

    cache_module.clear_all_view_caches()

    backend.clear.assert_called_once_with()
