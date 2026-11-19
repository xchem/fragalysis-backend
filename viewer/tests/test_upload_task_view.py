"""Tests for ``UploadTaskView`` reporting a finished computed-set upload.

A successful upload that collected warnings renders them as an HTML table via
pandas. pandas 2.x rejects ``display.max_colwidth=-1``, which used to turn that
response into a 500 even though the set had been loaded and committed.

The Celery result backend (``viewer.views.AsyncResult``) is patched, so no
broker/worker is needed.
"""
# pylint: disable=redefined-outer-name,unused-argument
import uuid
from unittest import mock

from django.test import Client
from django.urls import reverse

from viewer import views
from viewer.models import ComputedSet


def _mock_upload_result(monkeypatch, results):
    result = mock.Mock()
    result.status = "SUCCESS"
    result.id = str(uuid.uuid4())
    result.get.return_value = results
    monkeypatch.setattr(views, "AsyncResult", mock.Mock(return_value=result))


def test_successful_upload_with_warnings_returns_warning_table(
    monkeypatch, user, make_project, make_target
):
    project = make_project("proposal", members=[user])
    target = make_target(project)
    cset = ComputedSet.objects.create(name="cset", target=target, owner_user=user)

    process_messages = {
        "molecule_name": ["mol1"],
        "field": ["method"],
        "warning_string": ["Property method value docking ignored."],
    }
    _mock_upload_result(monkeypatch, ["process", cset.pk, process_messages])

    client = Client()
    client.force_login(user)
    url = reverse("viewer:upload_task", kwargs={"upload_task_id": uuid.uuid4()})
    response = client.get(url)

    assert response.status_code == 200, response.content
    data = response.json()
    assert data["validated"] == "Validated"
    assert "uploaded with warnings" in data["html"]
    assert "Property method value docking ignored." in data["html"]
