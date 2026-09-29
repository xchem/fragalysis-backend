"""Cross-origin isolation headers (issue #1041).

The frontend's Moorhen viewer uses threaded WebAssembly, which the browser only
permits when the page is *cross-origin isolated* (``window.crossOriginIsolated``).
That needs both of these headers on the top-level HTML document: -

    Cross-Origin-Opener-Policy: same-origin
    Cross-Origin-Embedder-Policy: require-corp

Django's ``SecurityMiddleware`` already provides COOP; the backend must add COEP.
Error responses are covered too, so a redirect or 404 doesn't lose isolation.
"""

import pytest
from webpack_loader.loader import WebpackLoader

_COOP = "Cross-Origin-Opener-Policy"
_COEP = "Cross-Origin-Embedder-Policy"


def _assert_isolated(response):
    assert response[_COOP] == "same-origin"
    assert response[_COEP] == "require-corp"


@pytest.fixture
def no_frontend_bundle(monkeypatch):
    """Render the React page without a built frontend (no webpack-stats.json)."""
    monkeypatch.setattr(WebpackLoader, "get_bundle", lambda self, name: [])


@pytest.mark.django_db
@pytest.mark.usefixtures("no_frontend_bundle")
def test_react_page_is_cross_origin_isolated(client):
    """The document the frontend is served from carries both policy headers."""
    response = client.get("/viewer/react/landing")

    assert response.status_code == 200
    _assert_isolated(response)


@pytest.mark.django_db
def test_not_found_response_is_cross_origin_isolated(client):
    """Error responses carry both policy headers."""
    response = client.get("/no-such-page/")

    assert response.status_code == 404
    _assert_isolated(response)


@pytest.mark.django_db
def test_redirect_response_is_cross_origin_isolated(client):
    """Redirects (e.g. the site root to the landing page) carry both headers."""
    response = client.get("/")

    assert response.status_code in (301, 302)
    _assert_isolated(response)
