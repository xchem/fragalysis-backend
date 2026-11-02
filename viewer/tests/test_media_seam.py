"""A guard for the media storage seam (issue 1038).

Media files can be on the filesystem or in S3, so all media access goes
through ``viewer.media_storage``. This test fails when other code uses
MEDIA_ROOT, a FieldFile's ``.path`` or ``default_storage``.

The code that still uses them is listed in ``_ALLOWED`` with its exact count.
Each task 3 PR (see #1038) moves code to the seam and lowers (or removes) its
count. A new use fails the test; so does a count that is now too high, so that
the list always shows the work that is left. When the list is empty, every
media access goes through the seam.

A line that matches but is not media access (i.e. a URL's ``.path``) can end
with the comment ``# media-seam: not media``.
"""
import re
from collections import Counter
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_PACKAGES = ("viewer", "api", "media_serve")
_SEAM = "viewer/media_storage.py"
_MARKER = "media-seam: not media"

_PATTERNS = (
    re.compile(r"\bMEDIA_ROOT\b"),
    # "<name>.path" that is not os.path, sys.path, request.path or self.path
    # (and not a method call such as "storage.path(name)").
    re.compile(
        r"(?<!\bos)(?<!\bsys)(?<!\brequest)(?<!\bself)(?<=[\w\)\]])\.path\b(?!\s*\()"
    ),
    re.compile(r"\bdefault_storage\b"),
)

# The remaining uses, as {file: count}. Lower these as the work is done.
_ALLOWED = {
    "viewer/download_restore.py": 2,
    "viewer/management/commands/cleanup_media.py": 3,
    "viewer/management/commands/restore_download_links.py": 2,
    "viewer/management/commands/restore_v3_files_from_v2_archive.py": 3,
    "viewer/management/commands/tags_from_sites.py": 1,
    "viewer/media_cleanup.py": 7,
    "viewer/media_watcher.py": 2,
    "viewer/models.py": 2,
    "viewer/squonk_job_file_transfer.py": 1,
    "viewer/target_delete.py": 1,
    "viewer/target_set_upload.py": 1,
    "viewer/upload_delete.py": 4,
    "viewer/utils.py": 4,
}


def _uses() -> Counter:
    uses: Counter = Counter()
    for package in _PACKAGES:
        for path in sorted((_ROOT / package).rglob("*.py")):
            name = path.relative_to(_ROOT).as_posix()
            if "/migrations/" in name or "/tests/" in name or name == _SEAM:
                continue
            for line in path.read_text(encoding="utf-8").splitlines():
                if _MARKER in line:
                    continue
                code = line.split("#", 1)[0]
                uses[name] += sum(len(p.findall(code)) for p in _PATTERNS)
    return +uses  # Removes the files with no uses


def test_media_access_goes_through_the_seam():
    uses = _uses()

    new = {f: n for f, n in uses.items() if n > _ALLOWED.get(f, 0)}
    assert not new, (
        f"New MEDIA_ROOT, FieldFile .path or default_storage uses: {new}. "
        "Use viewer.media_storage instead."
    )


def test_the_allowed_list_is_up_to_date():
    uses = _uses()

    lower = {f: uses.get(f, 0) for f, n in _ALLOWED.items() if uses.get(f, 0) < n}
    assert (
        not lower
    ), f"Media access moved to the seam. Lower (or remove) these in _ALLOWED: {lower}"
