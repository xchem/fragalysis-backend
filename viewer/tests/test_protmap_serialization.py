"""Regression test for ticket #958.

``SiteObservation.event_file`` holds a *binary* electron-density map. The
``ProtMapInfoSerializer`` formerly read it with
``open(..., encoding='utf-8').read()``, which raised ``UnicodeDecodeError`` on
the first non-UTF-8 byte (``0x9c`` in the reported traceback) and surfaced as a
500 on ``/api/protmap/``.

The serializer now returns the file's bytes base64-encoded, so the JSON
response is valid and lossless regardless of the file's contents.
"""
import base64
from pathlib import Path

from viewer.models import SiteObservation
from viewer.serializers import ProtMapInfoSerializer

# A map file is binary; this deliberately includes 0x9c (the byte from the
# ticket's traceback) which is not a valid UTF-8 start byte.
BINARY_MAP_CONTENT = b"MAP \x00\x9c\xff\xfe binary electron-density data"


def _write_event_file(media_root: str, content: bytes) -> str:
    """Write ``content`` under MEDIA_ROOT and return the storage-relative name."""
    name = "target_loader_data/test_958_event.map"
    target = Path(media_root) / name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(content)
    return name


def test_get_map_data_base64_encodes_binary_event_file(settings):
    """A binary event_file is returned base64-encoded, not decoded as text."""
    name = _write_event_file(settings.MEDIA_ROOT, BINARY_MAP_CONTENT)
    obs = SiteObservation(event_file=name)

    result = ProtMapInfoSerializer().get_map_data(obs)

    # It is valid base64 that round-trips back to the original bytes ...
    assert result == base64.b64encode(BINARY_MAP_CONTENT).decode("ascii")
    assert base64.b64decode(result) == BINARY_MAP_CONTENT
    # ... and crucially it is a str (JSON-renderable), never raising.
    assert isinstance(result, str)


def test_get_map_data_returns_none_when_no_event_file():
    """No event_file yields None rather than an error."""
    obs = SiteObservation(event_file="")

    assert ProtMapInfoSerializer().get_map_data(obs) is None
