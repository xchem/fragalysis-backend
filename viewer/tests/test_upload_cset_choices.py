"""Proposal/Target type-ahead suggestions on the ``viewer/upload_cset`` page.

The GET view offers, as ``target_choices``, only the targets whose proposal
the POST handler would accept from the user.
"""
# pylint: disable=redefined-outer-name,unused-argument
import pytest

URL = "/viewer/upload_cset/"


@pytest.fixture
def upload_settings(settings):
    """Force the authenticated-upload path with Django-only membership."""
    settings.AUTHENTICATE_UPLOAD = True
    settings.TA_AUTH_SERVICE = ""
    return settings


@pytest.fixture
def two_proposals(user, make_project, make_target):
    """One proposal the user is a member of, and one they are not."""
    mine = make_project("lb00001-1", members=[user])
    make_target(mine, title="Mpro")
    make_target(mine, title="Alpha")
    other = make_project("lb00002-1")
    make_target(other, title="Secret")


@pytest.mark.django_db
def test_choices_limited_to_users_proposals(
    upload_settings, two_proposals, authenticated_client
):
    response = authenticated_client.get(URL)

    assert response.status_code == 200
    assert response.context["target_choices"] == [
        {"project": "lb00001-1", "target": "Alpha"},
        {"project": "lb00001-1", "target": "Mpro"},
    ]


@pytest.mark.django_db
def test_choices_include_open_proposals(
    upload_settings, two_proposals, make_project, make_target, authenticated_client
):
    make_target(make_project("lb00003-1", open_to_public=True), title="Public")

    response = authenticated_client.get(URL)

    assert {"project": "lb00003-1", "target": "Public"} in response.context[
        "target_choices"
    ]


@pytest.mark.django_db
def test_all_targets_offered_without_upload_authentication(
    settings, two_proposals, authenticated_client
):
    settings.AUTHENTICATE_UPLOAD = False

    response = authenticated_client.get(URL)

    assert [c["target"] for c in response.context["target_choices"]] == [
        "Alpha",
        "Mpro",
        "Secret",
    ]


@pytest.mark.django_db
def test_page_wires_inputs_to_datalists(
    upload_settings, two_proposals, authenticated_client
):
    content = authenticated_client.get(URL).content.decode()

    assert 'list="proposal-options"' in content
    assert 'list="target-options"' in content
    assert '<datalist id="target-options">' in content
    assert '<script id="target-choices" type="application/json">' in content
    assert "Secret" not in content
