"""Representative end-to-end API access-control test.

This drives a real ISPyBSafeQuerySet-based viewset (``TargetView`` at
``/api/targets/``) through the DRF stack and asserts that results are filtered
to the proposals the requesting user can access. It is the template every
future viewset test should copy.

The test settings leave ``TA_AUTH_SERVICE`` unset, so proposal membership is
resolved from ``Project.user_id`` (no external service involved). Tests that
need the TA path can set ``settings.TA_AUTH_SERVICE`` and use the
``mock_target_access`` fixture instead.
"""


def _titles(response):
    """Target titles from a (LimitOffset-paginated) list response."""
    return {row["title"] for row in response.data["results"]}


def test_target_list_filtered_to_member_proposals(
    authenticated_client, user, make_project, make_target
):
    """An authenticated user sees only targets in proposals they belong to."""
    make_target(make_project("members-only", members=[user]), title="VisibleTarget")
    make_target(make_project("other-proposal"), title="HiddenTarget")

    response = authenticated_client.get("/api/targets/")

    assert response.status_code == 200
    assert _titles(response) == {"VisibleTarget"}


def test_target_list_includes_public_targets(
    authenticated_client, user, make_project, make_target
):
    """Open (public) proposals are visible even without membership."""
    make_target(make_project("members-only", members=[user]), title="VisibleTarget")
    make_target(
        make_project("public-proposal", open_to_public=True), title="PublicTarget"
    )
    make_target(make_project("private-proposal"), title="HiddenTarget")

    response = authenticated_client.get("/api/targets/")

    assert response.status_code == 200
    assert _titles(response) == {"VisibleTarget", "PublicTarget"}


def test_target_list_anonymous_sees_only_public(api_client, make_project, make_target):
    """An unauthenticated request sees only open (public) targets."""
    make_target(
        make_project("public-proposal", open_to_public=True), title="PublicTarget"
    )
    make_target(make_project("private-proposal"), title="HiddenTarget")

    response = api_client.get("/api/targets/")

    assert response.status_code == 200
    assert _titles(response) == {"PublicTarget"}
