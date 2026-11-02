"""The computed-set validate-only mode is gone (m2ms #2421).

Validate (``submit_choice=V``) had drifted from the upload code and gave
misleading answers. A failed upload rolls back and reports its error, so upload
is the validation: ``U`` hands the task params straight to
``process_compound_set`` and there is no ``validate_task`` endpoint any more.

``process_compound_set.delay`` is patched, so no broker/worker is needed.
"""
# pylint: disable=redefined-outer-name,unused-argument
import os
from unittest import mock

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client
from django.urls import NoReverseMatch, reverse

from viewer import tasks, views


@pytest.fixture
def upload_client(settings, user, make_project, make_target):
    """A logged-in client with a project/target it may upload to."""
    settings.AUTHENTICATE_UPLOAD = False
    project = make_project("proposal", members=[user])
    make_target(project, title="Mpro")
    client = Client()
    client.force_login(user)
    return client


@pytest.fixture
def mock_process(monkeypatch):
    task = mock.Mock()
    task.delay.return_value = mock.Mock(id="11111111-1111-1111-1111-111111111111")
    task.delay.return_value.status = "PENDING"
    monkeypatch.setattr(views, "process_compound_set", task)
    return task


def _post(client, choice):
    return client.post(
        reverse("viewer:upload_cset"),
        {
            "proposal_ref": "proposal",
            "target_name": "Mpro",
            "update_set": "None",
            "submit_choice": choice,
            "sdf_file": SimpleUploadedFile("cset.sdf", b"not really an sdf"),
        },
    )


def test_validate_task_url_is_gone():
    with pytest.raises(NoReverseMatch):
        reverse(
            "viewer:validate_task",
            kwargs={"validate_task_id": "11111111-1111-1111-1111-111111111111"},
        )


def test_validate_choice_starts_no_task(upload_client, mock_process):
    response = _post(upload_client, "V")

    assert response.status_code == 200
    mock_process.delay.assert_not_called()
    assert b"var taskUrl" not in response.content


def test_upload_choice_sends_params_straight_to_process(upload_client, mock_process):
    response = _post(upload_client, "U")

    assert response.status_code == 200
    mock_process.delay.assert_called_once()
    (task_params,) = mock_process.delay.call_args.args
    assert os.path.basename(task_params["sdf_file"]).startswith("cset")
    assert "update" in task_params
    assert b"/viewer/upload_task/" in response.content


@pytest.mark.parametrize("sdf_file", [None, "", "/no/such/file.sdf"])
def test_process_compound_set_rejects_missing_sdf(db, sdf_file):
    with pytest.raises(ValueError, match="SD file"):
        tasks.process_compound_set({"user_id": None, "sdf_file": sdf_file, "target": 1})
