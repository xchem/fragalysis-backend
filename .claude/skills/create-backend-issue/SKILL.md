---
name: create-backend-issue
description: Copy a project-board issue into the fragalysis-backend repo. Use when the user asks to create a backend issue from a project-board issue, e.g. "Create an issue for project issue 1234".
---

Copy an issue from the **m2ms project board**
(`https://github.com/orgs/m2ms/projects/2`) into this project's repository
(`xchem/fragalysis-backend`).

## When to use

Trigger phrases such as:

- "Create an issue for project issue 1234"
- "Make a backend issue from board issue 1234"

The number the user gives is the **project-board issue number**.

## What it does

Given a project-board issue number, the skill:

1. Reads the source issue from `m2ms/fragalysis-frontend` (the board's planning
   issues live there).
2. Verifies the issue carries the **`fragalysis-backend`** label. If it does
   not, it stops and reports — nothing is created.
3. Creates a new issue in `xchem/fragalysis-backend` with the **same title**
   and a **copy of the body** (plus a footer linking back to the source).
4. Adds a comment on the source issue with a link to the new backend issue.

It is **idempotent**: a hidden marker is written into the comment, so running
it twice for the same issue refuses to create a duplicate.

## How to run it

Run the helper script with the project-board issue number:

```bash
.claude/skills/create-backend-issue/create-backend-issue.sh <project-issue-number>
```

For example, for "Create an issue for project issue 1234":

```bash
.claude/skills/create-backend-issue/create-backend-issue.sh 1234
```

On success it prints the source URL and the new backend issue URL — relay both
to the user. On any error it exits non-zero and prints the reason to stderr;
report that reason rather than retrying blindly.

## Notes and assumptions

- Requires the `gh` CLI (authenticated, with `repo` scope) and `jq`. No GitHub
  *project* scope is needed because the issue is resolved by repository, not by
  querying the board API.
- The source repo is fixed to `m2ms/fragalysis-frontend`. This is safe because
  the `fragalysis-backend` label only exists in that repo, so any board issue
  carrying it lives there. If the board's planning issues ever move to another
  repo, update `SOURCE_REPO` in `create-backend-issue.sh`.
- The script never deletes or edits anything; it only creates one issue and one
  comment. If the label is missing, the issue cannot be read, or it was already
  imported, it aborts without side effects.
