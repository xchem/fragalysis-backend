---
name: release-staging
description: Release the staging branch to production. Opens a staging→production PR, waits for CI, merges it (keeping staging), then publishes a YYYY.MM.ITERATION GitHub release of production. Use when the user asks to "release staging", "do a release", or "release the staging branch".
---

Promote everything currently on **`staging`** to **`production`** and publish a
tagged GitHub release for `xchem/fragalysis-backend`.

This is an outward-facing, hard-to-reverse workflow: it merges to a protected
branch and publishes a public release that triggers the production image build.
Run the steps in order, **stop and report on any failure** (never swallow an
error), and pause for confirmation at the two checkpoints called out below.

## Key facts about this repo

- Remote/repo: `xchem/fragalysis-backend`.
- Release direction: PR with **`--base production --head staging`**.
- `production` is only ever updated *from* `staging`, so the PR normally has no
  conflicts. Conflicts only appear if someone committed directly to `production`
  (e.g. a hotfix). Handle them, don't assume they won't happen.
- **`staging` is a protected branch and must never be deleted.** Never pass
  `--delete-branch` to `gh pr merge`; never `git push --delete origin staging`.
- CI triggers (see `.github/workflows/`):
  - `build-staging.yaml` runs on every commit to `staging` — these are the
    checks that appear on the PR's head commit and are what you wait on.
  - `build-production.yaml` runs on **pushing a `YYYY.MM.N` tag**. It does *not*
    run on the merge to `production`. **Publishing the release (which creates the
    tag) is what starts the production build** — so the production build happens
    *after* the release, not before it.
- Tag scheme: `YYYY.MM.ITERATION`, ITERATION restarting at 1 each calendar month
  (latest example: `2026.06.1`). Use the helper script to compute it — don't do
  the arithmetic by hand.

## Steps

### 1. Get the Release title

Ask the user: **"What is the Release title for this release?"** Wait for their
answer. This exact string is used both as the PR title and the release title.
Do not invent one.

### 2. Pre-flight

```bash
git fetch origin staging production
```

Confirm there is actually something to release:

```bash
gh api repos/xchem/fragalysis-backend/compare/production...staging \
  --jq '{ahead: .ahead_by, behind: .behind_by, commits: [.commits[].commit.message]}'
```

If `ahead` is 0, there is nothing to release — stop and tell the user.

### 3. Create (or reuse) the PR

If an open `staging`→`production` PR already exists, reuse it; otherwise create
one with the title from step 1:

```bash
gh pr create --repo xchem/fragalysis-backend \
  --base production --head staging \
  --title "<Release title>" \
  --body "Release of staging to production."
```

Capture the PR number for the remaining steps. (To find an existing one:
`gh pr list --repo xchem/fragalysis-backend --base production --head staging --state open`.)

### 4. Resolve conflicts (only if any)

Check mergeability:

```bash
gh pr view <num> --repo xchem/fragalysis-backend --json mergeable,mergeStateStatus
```

If `mergeable` is `CONFLICTING`, resolve by bringing `production` into `staging`
(staging is the source of truth) — note `staging` is protected, so you may not be
able to push to it directly; if the push is rejected, **stop and ask the user how
they want conflicts resolved** rather than forcing anything. Re-check
mergeability before continuing.

### 5. Approval — normally not needed

`production` branch protection is configured with
`required_approving_review_count = 0`, so **no review approval is required** to
merge; do not try to approve the PR (GitHub forbids approving your *own* PR
anyway, so it would just error). Skip straight to waiting for CI.

Only if a merge is later blocked by a required-review rule (i.e. someone has
re-enabled required approvals) should you stop and ask the user to have a
*colleague* approve the PR — the author cannot self-approve. Confirm the current
rule with:

```bash
gh api repos/xchem/fragalysis-backend/branches/production/protection \
  --jq '.required_pull_request_reviews.required_approving_review_count'
```

### 6. Wait for a successful CI build

Watch the checks on the PR head and block until they finish:

```bash
gh pr checks <num> --repo xchem/fragalysis-backend --watch --fail-fast
```

If any check fails, **stop** and report which one — do not merge a red build.

### 7. Checkpoint, then merge (keep staging)

Once CI is green, summarise for the user (PR number/URL, that checks passed) and
**confirm before merging**. Then:

```bash
gh pr merge <num> --repo xchem/fragalysis-backend --merge
```

Use `--merge` (a real merge commit). **Never** add `--delete-branch` — `staging`
is protected and must survive. Afterwards verify it still exists:

```bash
git ls-remote --heads origin staging
```

### 8. Compute the release tag

```bash
.claude/skills/release-staging/next-release-tag.sh
```

It prints `NEXT_TAG=YYYY.MM.N` (the tag to create) and `PRIOR_TAG=...` (the
previous full release, used as the start point for the notes). It ignores
pre-release (`-rc.N`) tags.

### 9. Checkpoint, then publish the release

Show the user the computed `NEXT_TAG`, the title, and `PRIOR_TAG`, and **confirm
before publishing** (this is public and triggers the production build). Then:

```bash
gh release create <NEXT_TAG> --repo xchem/fragalysis-backend \
  --target production \
  --title "<Release title>" \
  --notes-start-tag <PRIOR_TAG> \
  --generate-notes
```

- `--target production` tags the current tip of `production` (the merge from
  step 7).
- `--generate-notes` with `--notes-start-tag <PRIOR_TAG>` makes the release notes
  contain **all changes since the prior release** (the merged PRs / commits
  between `PRIOR_TAG` and this tag) — exactly what's wanted. If `PRIOR_TAG` is
  empty (no previous release), omit `--notes-start-tag`.

### 10. Report

Tell the user, with links:

- the PR (and that it was merged, `staging` preserved),
- the new release URL and its tag,
- that pushing the tag has now started the `build production` CI, and they can
  watch it with
  `gh run list --repo xchem/fragalysis-backend --workflow "build production"`.

## Requirements

- `gh` CLI, authenticated with `repo` scope (and permission to merge to
  `production` / publish releases).
- `git`, `bash`, and the `gh`-bundled `jq` query engine (`--jq`).
