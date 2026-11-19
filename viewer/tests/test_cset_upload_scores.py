"""Tests for how computed-set uploads store score values (``MolOps``).

Each score column must end up as *one* ``ResultProperty`` per target, typed for
the whole column, rather than one per value type. The cases that used to split
a column:

- a missing value (``nan``), which was typed ``text``;
- an integer next to a float (``3`` and ``2.5``);
- integers at all, because ``dataType()`` returned ``"INT"`` while the stored
  type is called ``"integer"``, and ``0``/``1``, which were typed ``boolean``.

The uploads run through ``MolOps.task()`` against a hand-built target with one
observation that ``ref_mols`` and ``ref_pdb`` both resolve to.
"""
# pylint: disable=redefined-outer-name,unused-argument
import math
from datetime import datetime, timezone

import pytest
from rdkit import Chem

from viewer.cset_upload import MolOps, blank_mol_vals, column_data_type, value_data_type
from viewer.models import (
    Experiment,
    ExperimentUpload,
    Result,
    ResultProperty,
    ResultValueDataType,
    SiteObservation,
    TagCategory,
)

# the categories TagManager.tag_new_site_observations looks up
TAG_CATEGORIES = (
    "Other",
    "CanonSites",
    "ConformerSites",
    "Quatassemblies",
    "Crystalforms",
    "CrystalformSites",
)

REF_CODE = "Score-x0001_A_501_v1"

HEADER = {
    "ref_url": "https://example.com/method",
    "submitter_name": "Tester",
    "submitter_email": "tester@example.com",
    "submitter_institution": "Lab",
    "generation_date": "2026-10-27",
    "method": "scoretest",
}


@pytest.fixture
def score_target(settings, tmp_path, user, make_project, make_target):
    """A target with one upload and one observation to reference."""
    media_root = tmp_path / "media"
    media_root.mkdir()
    settings.MEDIA_ROOT = str(media_root)

    project = make_project("proposal", members=[user])
    target = make_target(project, title="Score")
    for category in TAG_CATEGORIES:
        TagCategory.objects.get_or_create(category=category)
    target.zip_archive = "Score_proposal"
    target.save()

    upload = ExperimentUpload.objects.create(
        project=project,
        target=target,
        committer=user,
        commit_datetime=datetime(2026, 1, 1, tzinfo=timezone.utc),
        upload_data_dir="upload_1",
        upload_version=1,
        data_version_major=3,
        data_version_minor=1,
        file="bundle_upload_1.tgz",
    )
    experiment = Experiment.objects.create(experiment_upload=upload, code="Score-x0001")
    SiteObservation.objects.create(
        code=REF_CODE, experiment=experiment, smiles="CCO", version=1
    )
    return target


@pytest.fixture
def upload_scores(tmp_path, user, score_target):
    """Upload an SDF whose data molecules carry the given score values.

    ``rows`` is a list of ``{score_name: value}`` dicts, one per molecule. The
    values are written as SDF text, as a submitter's file would hold them.
    """
    uploads = 0

    def _upload(rows):
        nonlocal uploads
        uploads += 1
        names = sorted({name for row in rows for name in row})

        sdf = tmp_path / f"upload_{uploads}.sdf"
        writer = Chem.SDWriter(str(sdf))
        header = Chem.MolFromSmiles("C")
        header.SetProp("_Name", "ver_1.2")
        for key, value in HEADER.items():
            header.SetProp(key, value)
        for name in names:
            header.SetProp(name, f"description of {name}")
        writer.write(header)

        for i, row in enumerate(rows):
            mol = Chem.MolFromSmiles("CC" + "C" * i)
            mol.SetProp("_Name", f"mol{uploads}_{i}")
            mol.SetProp("ref_mols", REF_CODE)
            mol.SetProp("ref_pdb", REF_CODE)
            for name, value in row.items():
                mol.SetProp(name, value)
            writer.write(mol)
        writer.close()

        submitter_name, submitter_method, version = blank_mol_vals(str(sdf))
        computed_set, _ = MolOps(
            user_id=user.pk,
            sdf_filename=str(sdf),
            submitter_name=submitter_name,
            submitter_method=submitter_method,
            target=score_target.pk,
            version=version,
            zfile=None,
            zfile_hashvals=None,
            computed_set_id=None,
        ).task()
        return computed_set

    return _upload


def _properties(target, name):
    return ResultProperty.objects.filter(target=target, result_property=name)


def _results(computed_set, name):
    return Result.objects.filter(
        computed_set=computed_set, result_property__result_property=name
    )


@pytest.mark.parametrize(
    "value, expected",
    [
        (3, "integer"),
        (0, "integer"),
        (1, "integer"),
        (2.5, "float"),
        (1.0, "float"),
        ("abc", "text"),
        ("True", "boolean"),
    ],
)
def test_value_data_type(value, expected):
    assert value_data_type(value) == expected


@pytest.mark.parametrize(
    "values, expected",
    [
        ([1, 2], "integer"),
        ([1, 2.5], "float"),
        ([1, "abc"], "text"),
        (["True", "False"], "boolean"),
        ([], None),
    ],
)
def test_column_data_type(values, expected):
    assert column_data_type(values) == expected


@pytest.mark.django_db
def test_missing_and_integer_values_keep_a_float_score_together(
    score_target, upload_scores
):
    """`nan`, an empty value and an int-looking value don't split a float score."""
    cset = upload_scores(
        [{"score": "1.5"}, {"score": "nan"}, {"score": "3"}, {"score": ""}]
    )

    properties = _properties(score_target, "score")
    assert properties.count() == 1
    assert properties.get().data_type.data_type == "float"

    results = _results(cset, "score")
    assert sorted(r.float_value for r in results) == [1.5, 3.0]
    assert all(r.text_value is None for r in results)


@pytest.mark.django_db
def test_integer_scores_are_stored_as_integers(score_target, upload_scores):
    """Including 0 and 1, which used to be taken for booleans."""
    cset = upload_scores([{"count": "1"}, {"count": "0"}, {"count": "7"}])

    assert _properties(score_target, "count").get().data_type.data_type == "integer"
    assert sorted(r.int_value for r in _results(cset, "count")) == [0, 1, 7]


@pytest.mark.django_db
def test_text_scores_are_stored_as_text(score_target, upload_scores):
    """One non-numeric value makes the whole column text, numbers included."""
    cset = upload_scores([{"label": "active"}, {"label": "2"}])

    assert _properties(score_target, "label").get().data_type.data_type == "text"
    assert sorted(r.text_value for r in _results(cset, "label")) == ["2", "active"]


@pytest.mark.django_db
def test_later_upload_reuses_the_existing_property(score_target, upload_scores):
    """An integer column is stored in an existing float property, not a new one."""
    upload_scores([{"score": "1.5"}])
    upload_scores([{"score": "2"}, {"score": "4"}])

    score = _properties(score_target, "score").get()
    assert score.data_type.data_type == "float"
    values = sorted(r.float_value for r in Result.objects.filter(result_property=score))
    assert values == [1.5, 2.0, 4.0]


@pytest.mark.django_db
def test_values_that_dont_fit_the_existing_property_fail_the_upload(
    score_target, upload_scores
):
    upload_scores([{"count": "1"}])

    with pytest.raises(ValueError, match="count"):
        upload_scores([{"count": "2.5"}])

    assert _properties(score_target, "count").count() == 1


@pytest.mark.django_db
def test_duplicate_properties_from_earlier_uploads_are_tolerated(
    score_target, upload_scores
):
    """Targets already split by the old bug: reuse the property of matching type."""
    for data_type in ("text", "float"):
        ResultProperty.objects.create(
            result_property="score",
            target=score_target,
            data_type=ResultValueDataType.objects.get(data_type=data_type),
        )

    cset = upload_scores([{"score": "1.5"}])

    assert _properties(score_target, "score").count() == 2
    (result,) = _results(cset, "score")
    assert result.result_property.data_type.data_type == "float"
    assert math.isclose(result.float_value, 1.5)


@pytest.mark.django_db
def test_duplicate_properties_prefer_one_the_values_fit(score_target, upload_scores):
    """An integer column goes into the float duplicate, not the older text one."""
    for data_type in ("text", "float"):
        ResultProperty.objects.create(
            result_property="score",
            target=score_target,
            data_type=ResultValueDataType.objects.get(data_type=data_type),
        )

    cset = upload_scores([{"score": "2"}])

    (result,) = _results(cset, "score")
    assert result.result_property.data_type.data_type == "float"
    assert result.float_value == 2.0
