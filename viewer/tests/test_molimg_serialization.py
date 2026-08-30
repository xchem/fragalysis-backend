"""Regression test for ticket #1011.

``MolImageSerializer.get_mol_image`` sized the image from the request's
``height``/``width`` query parameters, but guarded that lookup with a bare
``if params:`` - true whenever *any* parameter is present. A request that
filters or paginates without asking for a size, e.g. ``/api/molimg/?target=1``,
therefore raised ``MultiValueDictKeyError: 'height'`` and surfaced as a 500.

Each dimension now falls back to its default independently, so only the
supplied parameters influence the image.
"""
import re

import pytest
from rest_framework.request import Request
from rest_framework.test import APIRequestFactory

from viewer.models import SiteObservation
from viewer.serializers import MolImageSerializer

# Any valid SMILES will do; draw_mol() renders it to SVG.
SMILES = "c1ccccc1"

DEFAULT_HEIGHT = 125
DEFAULT_WIDTH = 125


def _mol_image(**query_params) -> str:
    """Serialize a SiteObservation's image for a request with these params."""
    request = Request(APIRequestFactory().get("/api/molimg/", query_params))
    serializer = MolImageSerializer(context={"request": request})
    return serializer.get_mol_image(SiteObservation(smiles=SMILES))


def _svg_dimensions(svg: str) -> tuple[int, int]:
    """Pull the (height, width) that RDKit wrote into the SVG element."""
    height = re.search(r"height=['\"](\d+)px['\"]", svg)
    width = re.search(r"width=['\"](\d+)px['\"]", svg)
    assert height and width, f"no dimensions in SVG: {svg[:200]}"
    return int(height.group(1)), int(width.group(1))


def test_get_mol_image_uses_defaults_when_size_params_absent():
    """A filtering param without height/width renders at the default size."""
    # Before the fix this raised MultiValueDictKeyError: 'height'.
    svg = _mol_image(target="1")

    assert _svg_dimensions(svg) == (DEFAULT_HEIGHT, DEFAULT_WIDTH)


def test_get_mol_image_uses_defaults_when_no_params_at_all():
    """The unparameterised request keeps rendering at the default size."""
    svg = _mol_image()

    assert _svg_dimensions(svg) == (DEFAULT_HEIGHT, DEFAULT_WIDTH)


def test_get_mol_image_honours_explicit_size_params():
    """Supplied height and width are used."""
    svg = _mol_image(height="200", width="300")

    assert _svg_dimensions(svg) == (200, 300)


@pytest.mark.parametrize(
    "params, expected",
    [
        ({"height": "200"}, (200, DEFAULT_WIDTH)),
        ({"width": "300"}, (DEFAULT_HEIGHT, 300)),
    ],
)
def test_get_mol_image_defaults_the_missing_dimension(params, expected):
    """One dimension given: the other falls back to its default."""
    assert _svg_dimensions(_mol_image(**params)) == expected


def test_get_mol_image_accepts_float_size_params():
    """Sizes arriving as floats are truncated to int, as they always were."""
    svg = _mol_image(height="200.7", width="300.2")

    assert _svg_dimensions(svg) == (200, 300)
