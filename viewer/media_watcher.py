"""Quieten the media deletion watcher while a deletion we asked for is running.

``filewatcher.sh`` logs every file removed under ``media/target_loader_data``. It
exists to catch files disappearing when nothing should have touched them, so the
deletions a user explicitly asked for - deleting a target, or peeling an upload
off it - are pure noise there, and a target deletion can bury the log in
thousands of lines.

Rather than stopping and restarting the watcher (which risks leaving it dead if
the deletion fails halfway), a deliberate deletion drops a marker file that the
watcher checks before writing each line. The marker carries the pid so
concurrent deletions cannot clear each other's, and it lives at MEDIA_ROOT -
outside the watched directory - so creating and removing it does not itself
generate events.

Suppression is best-effort on purpose: if the marker cannot be written the
deletion still goes ahead and the log simply gets its usual noise.
"""

import logging
import os
from contextlib import contextmanager
from pathlib import Path

from django.conf import settings

logger = logging.getLogger(__name__)

#: Marker files are ``<MEDIA_ROOT>/.filewatcher-pause.<pid>``; filewatcher.sh
#: looks for this prefix and stays quiet while any of them exists.
PAUSE_PREFIX = ".filewatcher-pause."


def _marker_path() -> Path:
    return Path(settings.MEDIA_ROOT) / f"{PAUSE_PREFIX}{os.getpid()}"


@contextmanager
def deletions_expected(reason: str):
    """Suppress the watcher's log for deletions made inside this block.

    :param reason: recorded in the marker, so a stale one left by a crashed
        process says what was running when it was written.
    """
    marker = _marker_path()
    created = False
    try:
        marker.write_text(f"{reason} (pid {os.getpid()})\n", encoding="utf-8")
        created = True
        logger.debug("Media deletion watcher paused: %s", reason)
    except OSError as exc:
        # Not fatal: the deletion matters, the log tidiness does not.
        logger.warning("Could not pause the media deletion watcher: %s", exc)

    try:
        yield
    finally:
        if created:
            try:
                marker.unlink()
                logger.debug("Media deletion watcher resumed")
            except OSError as exc:
                logger.warning(
                    "Could not resume the media deletion watcher; remove %s by "
                    "hand or it will stay quiet: %s",
                    marker,
                    exc,
                )
