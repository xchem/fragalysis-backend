"""A guard for the media storage seam (issue 1038).

Media files can be on the filesystem or in S3, so all media access goes
through ``viewer.media_storage``. This test fails when other code uses
MEDIA_ROOT, a FieldFile's ``.path`` or ``default_storage``.

The task 3 PRs (see #1038) moved all such code to the seam, so ``_ALLOWED`` is
empty. Keep it that way: a new use fails the test.

``_FILESYSTEM_ONLY`` lists the modules that are about the filesystem media
storage itself, with the reason. They are not scanned.

A line that matches but is not media access (i.e. a URL's ``.path``) can end
with the comment ``# media-seam: not media``.
"""
import re
from collections import Counter
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
_PACKAGES = ("viewer", "api", "media_serve")
_SEAM = "viewer/media_storage.py"

# Modules that only work with the filesystem media storage, by design.
_FILESYSTEM_ONLY = {
    # Quietens filewatcher.sh, which watches MEDIA_ROOT. It does nothing with S3
    # media storage, where the watcher is not started (see launch-stack.sh).
    "viewer/media_watcher.py",
}
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

# Allowed uses, as {file: count}. Empty: all media access goes through the seam.
_ALLOWED: dict[str, int] = {}


def _uses() -> Counter:
    uses: Counter = Counter()
    for package in _PACKAGES:
        for path in sorted((_ROOT / package).rglob("*.py")):
            name = path.relative_to(_ROOT).as_posix()
            if "/migrations/" in name or "/tests/" in name or name == _SEAM:
                continue
            if name in _FILESYSTEM_ONLY:
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


def test_the_filesystem_only_modules_exist():
    assert all((_ROOT / name).is_file() for name in _FILESYSTEM_ONLY)
