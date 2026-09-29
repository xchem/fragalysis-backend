"""Response headers nginx adds to static files (issue #1041).

Static files (``/static/``) are served directly by nginx, not Django, so the
Django cross-origin isolation middleware never sees them. Moorhen starts Web
Workers from static scripts (e.g. ``CootWorker.js``) and, because the page is
cross-origin isolated, each worker script needs its own
``Cross-Origin-Embedder-Policy: require-corp`` header or the browser blocks it.

An ``add_header`` in an nginx ``location`` *replaces* the ``add_header``s
inherited from the ``http`` block (it does not extend them), so the location
must also repeat every header declared in ``nginx.conf``.
"""

import re
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _location_block(conf: str, location: str) -> str:
    """The body of the named 'location' block in an nginx config."""
    match = re.search(r"location\s+" + re.escape(location) + r"\s*\{([^}]*)\}", conf)
    assert match, f"No 'location {location}' block found"
    return match.group(1)


def _added_headers(block: str) -> dict[str, str]:
    """Map of header name to 'add_header' line in an nginx block."""
    return {
        m.group(1): m.group(0)
        for m in re.finditer(r"^\s*add_header\s+(\S+)\s+[^;]*;", block, re.MULTILINE)
    }


def test_static_files_get_cross_origin_embedder_policy():
    conf = (_REPO_ROOT / "django_nginx.conf").read_text()
    static_headers = _added_headers(_location_block(conf, "/static/"))

    assert "Cross-Origin-Embedder-Policy" in static_headers
    assert '"require-corp"' in static_headers["Cross-Origin-Embedder-Policy"]


def test_static_files_keep_the_http_level_headers():
    """Declaring a header in the location must not drop the nginx.conf ones."""
    http_headers = _added_headers((_REPO_ROOT / "nginx.conf").read_text())
    conf = (_REPO_ROOT / "django_nginx.conf").read_text()
    static_headers = _added_headers(_location_block(conf, "/static/"))

    assert http_headers, "Expected add_header directives in nginx.conf"
    assert set(http_headers) <= set(static_headers)
