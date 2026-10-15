"""The ``/api/molimg/`` endpoint has been removed (issue 1011).

It rendered a 2D SVG for each SiteObservation, was not used by the front-end,
and failed with a 500 whenever a query parameter other than ``height``/``width``
was supplied. Rather than fix it, it was retired.
"""

import pytest


@pytest.mark.django_db
def test_molimg_endpoint_is_not_found(authenticated_client):
    """Requests to the old endpoint get a 404, not a 500."""
    response = authenticated_client.get("/api/molimg/", {"target": 1})

    assert response.status_code == 404


@pytest.mark.django_db
def test_molimg_is_not_listed_in_the_api_root(authenticated_client):
    """The browsable API root no longer advertises the endpoint."""
    response = authenticated_client.get("/api/")

    assert response.status_code == 200
    assert "molimg" not in response.data
