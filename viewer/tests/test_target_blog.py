"""Tests for the manually entered target blog card (m2ms-2228).

A ``TargetBlog`` holds the blog card shown next to a target: title, publishing
date, overview, image and the post's URL. It is written through
``/api/target-blogs/`` by members of the target's proposal and read by anyone
who can see the target. The image is stored in the target's own
``target_loader_data/`` subdirectory, so deleting the target removes it too.
It is served by its own media route, because the generic ``target_loader_data``
route sends a non-image content type.

``django_cleanup`` deletes files in ``transaction.on_commit``, and the ordinary
``db`` fixture never commits, so the file-removal tests use
``django_capture_on_commit_callbacks(execute=True)``.
"""

# Fixtures legitimately reuse their names as arguments, and `media_root` is
# requested for its side effect only - standard pytest patterns
# that pylint misreads.
# pylint: disable=redefined-outer-name,unused-argument
from io import BytesIO
from pathlib import Path

import pytest
import requests
from django.core.files.uploadedfile import SimpleUploadedFile
from PIL import Image

from viewer import target_blog
from viewer.media_cleanup import scan_target_loader_data
from viewer.models import TargetBlog, UserRole
from viewer.target_delete import delete_target

_ENDPOINT = "/api/target-blogs/"


def _png(
    name: str = "blog.png", colour: str = "red", size: tuple[int, int] = (4, 4)
) -> SimpleUploadedFile:
    buffer = BytesIO()
    Image.new("RGB", size, colour).save(buffer, format="PNG")
    return SimpleUploadedFile(name, buffer.getvalue(), content_type="image/png")


@pytest.fixture
def media_root(settings, tmp_path) -> Path:
    root = tmp_path / "media"
    (root / settings.TARGET_LOADER_MEDIA_DIRECTORY).mkdir(parents=True)
    settings.MEDIA_ROOT = str(root)
    return root


@pytest.fixture
def member(make_user):
    return make_user("member")


@pytest.fixture
def outsider(make_user):
    return make_user("outsider")


@pytest.fixture
def make_loaded_target(make_project, make_target, member, settings):
    """A Target that has been through the loader, so it owns a media subdir."""
    settings.TA_AUTH_SERVICE = ""

    def _make(title="Mpro", proposal="lb12345-1", open_to_public=False):
        project = make_project(
            proposal, members=[member], open_to_public=open_to_public
        )
        target = make_target(project, title=title)
        target.zip_archive = f"{title}_{proposal}"
        target.save()
        return target

    return _make


def _payload(target, **overrides):
    data = {
        "target": target.pk,
        "title": "Blog: Affinity and kinetics data",
        "published": "2026-05-05",
        "overview": "How this first data package was generated.",
        "url": "https://openbind.uk/news/blog-affinity-and-kinetics/",
    }
    data.update(overrides)
    return data


def _create(client, target, **overrides):
    return client.post(_ENDPOINT, _payload(target, **overrides), format="multipart")


# ---------------------------------------------------------------------------
# create


def test_member_creates_blog_with_image(
    api_client, member, make_loaded_target, media_root, settings
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image=_png())

    assert response.status_code == 201, response.data
    blog = TargetBlog.objects.get(target=target)
    assert blog.title == "Blog: Affinity and kinetics data"
    assert str(blog.published) == "2026-05-05"
    assert blog.image.name.startswith(
        f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/{target.zip_archive}/blog/"
    )
    assert (media_root / blog.image.name).is_file()
    assert response.data["image"].endswith(f"/media/target_blog/{blog.pk}/")


def test_image_and_optional_fields_can_be_omitted(
    api_client, member, make_loaded_target, media_root
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = api_client.post(
        _ENDPOINT,
        {"target": target.pk, "title": "A post", "url": "https://openbind.uk/n/a/"},
        format="multipart",
    )

    assert response.status_code == 201, response.data
    assert response.data["image"] is None
    assert response.data["published"] is None
    assert response.data["overview"] == ""


@pytest.mark.parametrize("missing", ["title", "url"])
def test_title_and_url_are_required(
    api_client, member, make_loaded_target, media_root, missing
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    data = _payload(target)
    del data[missing]

    response = api_client.post(_ENDPOINT, data, format="multipart")

    assert response.status_code == 400
    assert missing in response.data


def test_non_member_cannot_create(api_client, outsider, make_loaded_target, media_root):
    target = make_loaded_target()
    api_client.force_authenticate(user=outsider)

    response = _create(api_client, target)

    assert response.status_code in (400, 403)
    assert not TargetBlog.objects.exists()


def test_non_member_cannot_create_on_public_target(
    api_client, outsider, make_loaded_target, media_root
):
    """Seeing a public target is not the same as being allowed to edit it."""
    target = make_loaded_target(open_to_public=True)
    api_client.force_authenticate(user=outsider)

    response = _create(api_client, target)

    assert response.status_code in (400, 403)
    assert not TargetBlog.objects.exists()


def test_anonymous_cannot_create(api_client, make_loaded_target, media_root):
    target = make_loaded_target(open_to_public=True)

    response = _create(api_client, target)

    assert response.status_code in (400, 401, 403)
    assert not TargetBlog.objects.exists()


def test_second_blog_for_same_target_is_rejected(
    api_client, member, make_loaded_target, media_root
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    assert _create(api_client, target).status_code == 201

    response = _create(api_client, target, title="Another post")

    assert response.status_code == 400
    assert TargetBlog.objects.filter(target=target).count() == 1


def test_target_without_media_directory_rejects_image(
    api_client, member, make_loaded_target, media_root
):
    target = make_loaded_target()
    target.zip_archive = None
    target.save()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image=_png())

    assert response.status_code == 400
    assert not TargetBlog.objects.exists()
    assert not any(media_root.rglob("*.png"))


def test_non_image_file_is_rejected(api_client, member, make_loaded_target, media_root):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    fake = SimpleUploadedFile("x.png", b"not an image", content_type="image/png")

    response = _create(api_client, target, image=fake)

    assert response.status_code == 400
    assert "image" in response.data


def test_image_over_byte_limit_is_rejected(
    api_client, member, make_loaded_target, media_root, settings
):
    settings.TARGET_BLOG_IMAGE_MAX_BYTES = 10
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image=_png())

    assert response.status_code == 400
    assert "image" in response.data


@pytest.mark.parametrize("size", [(1500, 1000), (3000, 3000), (1000, 3000)])
def test_image_within_pixel_limit_is_accepted(
    api_client, member, make_loaded_target, media_root, size
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image=_png(size=size))

    assert response.status_code == 201, response.data


@pytest.mark.parametrize("size", [(3001, 1000), (1000, 3001)])
def test_image_over_pixel_limit_is_rejected(
    api_client, member, make_loaded_target, media_root, size
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image=_png(size=size))

    assert response.status_code == 400
    assert "image" in response.data
    assert not TargetBlog.objects.exists()


def test_without_upload_authentication_anyone_can_write(
    api_client, outsider, make_loaded_target, media_root, settings
):
    """AUTHENTICATE_UPLOAD=False (development only) drops the membership
    check, as it does for the other uploads and for target deletion.
    """
    settings.AUTHENTICATE_UPLOAD = False
    target = make_loaded_target(open_to_public=True)

    created = _create(api_client, target)
    assert created.status_code == 201, created.data
    blog_url = f"{_ENDPOINT}{created.data['id']}/"

    api_client.force_authenticate(user=outsider)
    patched = api_client.patch(blog_url, {"title": "New"}, format="multipart")
    assert patched.status_code == 200, patched.data
    assert api_client.delete(blog_url).status_code == 204


def test_loader_role_does_not_grant_blog_writes(
    api_client, outsider, make_loaded_target, media_root
):
    """The Loader role is for automated uploads; it skips the membership check
    there, but not for blog cards.
    """
    role, _ = UserRole.objects.get_or_create(name=UserRole.LOADER_ROLE)
    role.users.add(outsider)
    target = make_loaded_target(open_to_public=True)
    blog = TargetBlog.objects.create(target=target, title="Old", url="https://x/p/")
    other = make_loaded_target(title="Other", proposal="lb12345-2", open_to_public=True)
    api_client.force_authenticate(user=outsider)

    assert _create(api_client, other).status_code == 403
    assert not TargetBlog.objects.filter(target=other).exists()
    patched = api_client.patch(
        f"{_ENDPOINT}{blog.pk}/", {"title": "New"}, format="multipart"
    )
    assert patched.status_code == 403
    assert api_client.delete(f"{_ENDPOINT}{blog.pk}/").status_code == 403
    blog.refresh_from_db()
    assert blog.title == "Old"


# ---------------------------------------------------------------------------
# create from image_url

_IMAGE_URL = "https://openbind.uk/wp-content/uploads/2026/04/EV2A.png"


class _FakeResponse:
    def __init__(self, body=b"", status_code=200, is_redirect=False):
        self.body = body
        self.status_code = status_code
        self.is_redirect = is_redirect
        self.headers = {"Content-Type": "image/png"}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def iter_content(self, chunk_size):
        for start in range(0, len(self.body), chunk_size):
            yield self.body[start : start + chunk_size]


class _FakeGet:
    """Stands in for requests.get. Records the fetched URLs; set ``response``
    (or ``error``) to control the answer.
    """

    def __init__(self):
        self.calls: list = []
        self.response = _FakeResponse(_png().read())
        self.error: Exception | None = None

    def __call__(self, url, **kwargs):
        self.calls.append((url, kwargs))
        if self.error is not None:
            raise self.error
        return self.response


@pytest.fixture
def fake_get(monkeypatch) -> _FakeGet:
    fake = _FakeGet()
    monkeypatch.setattr(target_blog.requests, "get", fake)
    return fake


def test_member_creates_blog_from_image_url(
    api_client, member, make_loaded_target, media_root, fake_get, settings
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url=_IMAGE_URL)

    assert response.status_code == 201, response.data
    assert "image_url" not in response.data
    blog = TargetBlog.objects.get(target=target)
    assert blog.image.name == (
        f"{settings.TARGET_LOADER_MEDIA_DIRECTORY}/{target.zip_archive}/blog/EV2A.png"
    )
    assert (media_root / blog.image.name).is_file()
    [(url, kwargs)] = fake_get.calls
    assert url == _IMAGE_URL
    assert kwargs["allow_redirects"] is False
    assert kwargs["timeout"]


def test_image_url_without_extension_gets_one(
    api_client, member, make_loaded_target, media_root, fake_get
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url="https://openbind.uk/img/42")

    assert response.status_code == 201, response.data
    assert TargetBlog.objects.get(target=target).image.name.endswith("/blog/42.png")


def test_member_replaces_image_from_url(
    api_client, member, make_loaded_target, media_root, fake_get
):
    target = make_loaded_target()
    blog = TargetBlog.objects.create(target=target, title="Post", url="https://x/p/")
    api_client.force_authenticate(user=member)

    response = api_client.patch(
        f"{_ENDPOINT}{blog.pk}/", {"image_url": _IMAGE_URL}, format="multipart"
    )

    assert response.status_code == 200, response.data
    blog.refresh_from_db()
    assert blog.image.name.endswith("/blog/EV2A.png")


@pytest.mark.parametrize(
    "url",
    [
        "https://example.com/a.png",
        "http://169.254.169.254/latest/meta-data/",
        "https://openbind.uk.evil.example/a.png",
        "ftp://openbind.uk/a.png",
    ],
)
def test_image_url_off_allow_list_is_never_fetched(
    api_client, member, make_loaded_target, media_root, fake_get, url
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url=url)

    assert response.status_code == 400
    assert "image_url" in response.data
    assert fake_get.calls == []
    assert not TargetBlog.objects.exists()


def test_empty_allow_list_disables_image_url(
    api_client, member, make_loaded_target, media_root, fake_get, settings
):
    settings.TARGET_BLOG_IMAGE_URL_ALLOWED_HOSTS_LIST = []
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url=_IMAGE_URL)

    assert response.status_code == 400
    assert fake_get.calls == []


def test_non_member_cannot_make_backend_fetch(
    api_client, outsider, make_loaded_target, media_root, fake_get
):
    target = make_loaded_target(open_to_public=True)
    api_client.force_authenticate(user=outsider)

    response = _create(api_client, target, image_url=_IMAGE_URL)

    assert response.status_code in (400, 403)
    assert fake_get.calls == []


def test_image_and_image_url_together_are_rejected(
    api_client, member, make_loaded_target, media_root, fake_get
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image=_png(), image_url=_IMAGE_URL)

    assert response.status_code == 400
    assert "image_url" in response.data
    assert fake_get.calls == []


@pytest.mark.parametrize(
    "answer",
    [
        _FakeResponse(status_code=404),
        _FakeResponse(status_code=301, is_redirect=True),
        _FakeResponse(b"<html>not an image</html>"),
        requests.exceptions.ConnectTimeout("too slow"),
    ],
    ids=["404", "redirect", "not-an-image", "timeout"],
)
def test_unusable_image_url_is_a_400(
    api_client, member, make_loaded_target, media_root, fake_get, answer
):
    if isinstance(answer, Exception):
        fake_get.error = answer
    else:
        fake_get.response = answer
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url=_IMAGE_URL)

    assert response.status_code == 400
    assert "image_url" in response.data
    assert not TargetBlog.objects.exists()


def test_image_url_download_stops_at_byte_limit(
    api_client, member, make_loaded_target, media_root, fake_get, settings
):
    settings.TARGET_BLOG_IMAGE_MAX_BYTES = 100
    fake_get.response = _FakeResponse(b"x" * 10_000_000)
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url=_IMAGE_URL)

    assert response.status_code == 400
    assert "larger than 100 bytes" in str(response.data["image_url"])


def test_image_url_over_pixel_limit_is_rejected(
    api_client, member, make_loaded_target, media_root, fake_get
):
    fake_get.response = _FakeResponse(_png(size=(3001, 10)).read())
    target = make_loaded_target()
    api_client.force_authenticate(user=member)

    response = _create(api_client, target, image_url=_IMAGE_URL)

    assert response.status_code == 400
    assert "image_url" in response.data


# ---------------------------------------------------------------------------
# read


def test_anonymous_reads_public_targets_blog(
    api_client, make_loaded_target, media_root
):
    target = make_loaded_target(open_to_public=True)
    TargetBlog.objects.create(target=target, title="Post", url="https://x.org/p/")

    response = api_client.get(_ENDPOINT, {"target": target.pk})

    assert response.status_code == 200
    assert [row["title"] for row in response.data["results"]] == ["Post"]


def test_non_member_does_not_see_private_targets_blog(
    api_client, outsider, make_loaded_target, media_root
):
    target = make_loaded_target()
    blog = TargetBlog.objects.create(target=target, title="Post", url="https://x/p/")
    api_client.force_authenticate(user=outsider)

    assert api_client.get(_ENDPOINT).data["results"] == []
    assert api_client.get(f"{_ENDPOINT}{blog.pk}/").status_code == 404


# ---------------------------------------------------------------------------
# update / delete


def test_member_updates_fields(api_client, member, make_loaded_target, media_root):
    target = make_loaded_target()
    blog = TargetBlog.objects.create(target=target, title="Old", url="https://x/p/")
    api_client.force_authenticate(user=member)

    response = api_client.patch(
        f"{_ENDPOINT}{blog.pk}/", {"title": "New"}, format="multipart"
    )

    assert response.status_code == 200, response.data
    blog.refresh_from_db()
    assert blog.title == "New"


def test_target_cannot_be_changed_on_update(
    api_client, member, make_loaded_target, media_root
):
    target = make_loaded_target()
    other = make_loaded_target(title="Other", proposal="lb12345-2")
    blog = TargetBlog.objects.create(target=target, title="Post", url="https://x/p/")
    api_client.force_authenticate(user=member)

    api_client.patch(f"{_ENDPOINT}{blog.pk}/", {"target": other.pk}, format="multipart")

    blog.refresh_from_db()
    assert blog.target == target


def test_non_member_cannot_update_public_targets_blog(
    api_client, outsider, make_loaded_target, media_root
):
    target = make_loaded_target(open_to_public=True)
    blog = TargetBlog.objects.create(target=target, title="Old", url="https://x/p/")
    api_client.force_authenticate(user=outsider)

    response = api_client.patch(
        f"{_ENDPOINT}{blog.pk}/", {"title": "New"}, format="multipart"
    )

    assert response.status_code == 403
    blog.refresh_from_db()
    assert blog.title == "Old"


def test_replacing_image_removes_old_file(
    api_client,
    member,
    make_loaded_target,
    media_root,
    django_capture_on_commit_callbacks,
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    with django_capture_on_commit_callbacks(execute=True):
        created = _create(api_client, target, image=_png("one.png"))
    old_file = media_root / TargetBlog.objects.get(pk=created.data["id"]).image.name
    assert old_file.is_file()

    with django_capture_on_commit_callbacks(execute=True):
        response = api_client.patch(
            f"{_ENDPOINT}{created.data['id']}/",
            {"image": _png("two.png", colour="blue")},
            format="multipart",
        )

    assert response.status_code == 200, response.data
    new_file = media_root / TargetBlog.objects.get(pk=created.data["id"]).image.name
    assert new_file.is_file()
    assert not old_file.exists()


def test_member_deletes_blog_and_its_file(
    api_client,
    member,
    make_loaded_target,
    media_root,
    django_capture_on_commit_callbacks,
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    with django_capture_on_commit_callbacks(execute=True):
        created = _create(api_client, target, image=_png())
    image = media_root / TargetBlog.objects.get(pk=created.data["id"]).image.name

    with django_capture_on_commit_callbacks(execute=True):
        response = api_client.delete(f"{_ENDPOINT}{created.data['id']}/")

    assert response.status_code == 204
    assert not TargetBlog.objects.exists()
    assert not image.exists()


def test_delete_target_removes_blog_and_image(
    api_client,
    member,
    make_loaded_target,
    media_root,
    django_capture_on_commit_callbacks,
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    with django_capture_on_commit_callbacks(execute=True):
        _create(api_client, target, image=_png())
    image = media_root / TargetBlog.objects.get(target=target).image.name

    with django_capture_on_commit_callbacks(execute=True):
        delete_target(target)

    assert not TargetBlog.objects.exists()
    assert not image.exists()


def test_cleanup_media_does_not_treat_blog_dir_as_orphaned(
    api_client, member, make_loaded_target, media_root
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    _create(api_client, target, image=_png())

    report = scan_target_loader_data()

    assert report.candidates == []


# ---------------------------------------------------------------------------
# image media route


def test_image_route_redirects_to_nginx_as_inline_image(
    api_client, member, make_loaded_target, media_root
):
    target = make_loaded_target(open_to_public=True)
    api_client.force_authenticate(user=member)
    blog_id = _create(api_client, target, image=_png()).data["id"]
    blog = TargetBlog.objects.get(pk=blog_id)
    api_client.force_authenticate(user=None)

    response = api_client.get(f"/media/target_blog/{blog_id}/")

    assert response.status_code == 200
    assert response["X-Accel-Redirect"] == f"/{blog.image.name}"
    assert response["Content-Type"] == "image/png"
    assert response["Content-Disposition"].startswith("inline")


def test_image_route_hides_private_targets_blog(
    api_client, member, outsider, make_loaded_target, media_root
):
    target = make_loaded_target()
    api_client.force_authenticate(user=member)
    blog_id = _create(api_client, target, image=_png()).data["id"]
    api_client.force_authenticate(user=None)

    # The media route is a plain Django view, so it needs a session login;
    # force_authenticate only reaches DRF views.
    api_client.force_login(member)
    assert api_client.get(f"/media/target_blog/{blog_id}/").status_code == 200
    api_client.force_login(outsider)
    assert api_client.get(f"/media/target_blog/{blog_id}/").status_code == 404


def test_image_route_404s_without_image(api_client, make_loaded_target, media_root):
    target = make_loaded_target(open_to_public=True)
    blog = TargetBlog.objects.create(target=target, title="Post", url="https://x/p/")

    assert api_client.get(f"/media/target_blog/{blog.pk}/").status_code == 404
