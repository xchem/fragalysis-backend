---
name: issue
description: Work on a GitHub issue end-to-end in a git worktree — move its board card to In Progress first, read the issue, explore the code, branch from staging, write the test first (TDD), commit, and open a PR. Use when the user says e.g. "work on issue 958", "start issue 1234", "/issue N", or asks to take a fragalysis-backend issue from triage to PR.
---

# Issue Workflow

Take a GitHub issue from triage to an open PR, working in an isolated git
**worktree**. The issue number is given as the argument (e.g. `958`).

## Instructions

### 1. Move the board card to "In Progress"

Do this **first**, before reading the issue (which can take some time).

The m2ms project board (`https://github.com/orgs/m2ms/projects/2`) tracks the
*frontend* source issue, not this backend issue. Move its card from its backlog
lane into the matching "In Progress" lane — e.g. `Infra - Backlog` →
`Infra - In Progress`:

```bash
.claude/skills/issue/move-board-card.sh <number> --in-progress
```

The script follows the backend issue's "source issue" footer to the board card,
derives the lane's category (the text before the first ` - `), and moves the card
to `{Category} - In Progress`. Cards **already in an "In Progress" lane are left
untouched.** If there is no matching "In Progress" lane — or no board card can be
found — it changes nothing and prints why; **relay that to the user**. Requires a
`gh` token with the `project` (write) scope.

### 2. Read the issue

Fetch and understand the issue from `xchem/fragalysis-backend`:

```bash
gh issue view <number> --repo xchem/fragalysis-backend \
  --json number,title,body,labels,state,comments
```

Summarise for the user: what needs doing, the acceptance criteria, and any
blockers or dependencies. Follow links to source/board issues if the body has
them (backend issues imported by `create-backend-issue` carry a footer linking
to their `m2ms/fragalysis-frontend` source).

### 3. Explore the codebase

Before coding: search for the related code, understand the current
implementation, and identify the files that need changing. Read
`ARCHITECTURE.md` first if the change touches models, security, or the loader.

### 4. Create a worktree

We always branch from **`staging`** (the default branch). The branch name is the
committer's initials, the issue number, and a brief description:

    {initials}-{issue-number}-{brief-description}

Create a worktree whose git branch is prefixed with the repository name, in a
sibling directory named for the branch:

```bash
git fetch origin
repo_name="$(basename "$(git remote get-url origin)" .git)"   # fragalysis-backend
branch="{initials}-{issue-number}-{brief-description}"
git worktree add -b "${repo_name}_${branch}" "../${branch}" staging
cd "../${branch}"
```

> **Poetry gotcha:** a new worktree is a new path, so Poetry creates a *fresh,
> empty* virtualenv for it — `poetry run pytest` there fails with
> `ModuleNotFoundError: No module named 'django'`. Either run
> `poetry install --only main,test` in the worktree, or reuse the main
> checkout's venv: get its path once with
> `poetry env info --path` (run from the main checkout) and invoke that
> `…/bin/pytest` directly from inside the worktree.

### 5. Implement the changes (TDD)

- **Write the test first**, watch it fail (red), then make it pass (green). See
  `viewer/tests/` for the fixture style (`conftest.py`), and existing regression
  tests like `test_experiment_serialization.py` for the pattern of a
  "500 → fixed" serializer bug.
- Follow project patterns and the conventions in `CLAUDE.md` (access control via
  `ISPyBSafeQuerySet`, environment-driven settings, never let errors pass
  silently).
- Run the suite against a real Postgres (SQLite can't stand in):
  `docker compose -f docker-compose.test.yml up -d --wait database` then pytest.
- Run `pre-commit run --files <changed files>` before committing (lint is limited
  to the `viewer` app).
- Make incremental, **Conventional Commit** messages that name the issue, e.g.
  `fix(958): base64-encode binary map data in /api/protmap/`. End commit messages
  with the trailer:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

### 6. Update the issue

As you work, comment on the issue with progress and any blockers/questions, and
note the PR once it exists.

### 7. Create the PR and link it

```bash
gh pr create --repo xchem/fragalysis-backend --base staging \
  --head "${repo_name}_${branch}" \
  --title "fix(<issue>): <description>" \
  --body-file <file>
```

- Put the issue number in the PR title (`feat(123): …` / `fix(123): …`).
- Link the issue: open the body with `Fixes #<number>` so it auto-closes on
  merge.
- **Write the PR body to a file and use `--body-file`.** Passing a long
  `--body "$(cat <<EOF …)"` lets the shell evaluate backticks/`$()` in the body
  and mangles it — a file avoids that.

Once the PR is open, move the board card into the review lane:

```bash
.claude/skills/issue/move-board-card.sh <number> "Dev Done - Do review (DEV)"
```

If that lane can't be found or no board card exists, the script changes nothing
and says why — **relay that to the user**.

### 8. If you find an unrelated bug

1. Create a new issue with the details (`create-backend-issue` if it originates
   on the board).
2. Link it to the current issue if related.
3. Note it in the PR description.
4. Continue with the original task — don't scope-creep.

## Notes

- The worktree is left in place for review; remove it after the PR merges with
  `git worktree remove ../{branch}`.
- `staging` is a protected branch — never push to it directly or delete it.
