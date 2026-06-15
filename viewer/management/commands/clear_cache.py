"""Drop every entry in the configured page-cache alias.

Useful after schema migrations, manual DB edits, or any out-of-band change
that wouldn't have triggered the per-view invalidation hooks in views.py
and tasks.py.
"""
from django.conf import settings
from django.core.cache import caches
from django.core.management.base import BaseCommand

from viewer.cache import clear_all_view_caches


class Command(BaseCommand):
    help = "Clear the entire page cache (configured CACHE_MIDDLEWARE_ALIAS)."

    def handle(self, *args, **kwargs):
        del args, kwargs
        alias = settings.CACHE_MIDDLEWARE_ALIAS
        backend = type(caches[alias]).__name__
        clear_all_view_caches()
        self.stdout.write(
            self.style.SUCCESS(f"Cleared cache alias '{alias}' ({backend}).")
        )
