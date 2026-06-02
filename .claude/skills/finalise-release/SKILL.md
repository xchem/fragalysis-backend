---
name: finalise-release
description: Finalise a Release ticket on the m2ms Fragalysis project board — move every "fragalysis-backend"-labelled issue from the "Approved in staging - push to production" swimlane into "In production (Done)" and list them on the Release ticket. Use after a staging→production release, when the user says e.g. "finalise release 2246", "mark release 2246 done", "move the approved backend issues to production for release 2246".
---

Finalise a **Release** ticket on the org-level Fragalysis project board
(`https://github.com/orgs/m2ms/projects/2`).

This is the board-bookkeeping step that normally follows the
[`release-staging`](../release-staging/SKILL.md) skill: once staging has been
pushed to production, the backend issues that were sitting in the
**"Approved in staging - push to production"** swimlane need to be marked
**"In production (Done)"**, and recorded on the Release ticket.

## When to use

The user supplies a **Release ticket's issue number** (a `m2ms/fragalysis-frontend`
issue number — that's where the board's planning/Release tickets live). Trigger
phrases such as:

- "Finalise release 2246"
- "Mark release 2246 as done / in production"
- "Move the approved backend issues to production for release 2246"

## What it does

Given the Release ticket's issue number, the skill:

1. Resolves the board's **Status** single-select field and **verifies the three
   required swimlanes exist before doing anything**:
   - `Releases`
   - `Approved in staging - push to production`
   - `In production (Done)`

   If any is missing it aborts without making changes.
2. **Validates the Release ticket**: the given issue must be on the board, sit in
   the `Releases` swimlane, and have a title starting with the word "Release".
   Otherwise it stops — this guards against finalising the wrong ticket.
3. Finds every issue in `Approved in staging - push to production` that carries
   the **`fragalysis-backend`** label.
4. **Moves** each of those issues to `In production (Done)`.
5. **Appends** a bullet list of the moved issues (issue IDs + titles, as links)
   to the Release ticket's description, under a
   "Backend issues released to production" heading.

If there is **nothing to move** it changes nothing and says so prominently —
relay that to the user, because it often means a mistake (wrong Release ticket,
or the issues were never approved into the swimlane).

## How to run it

```bash
.claude/skills/finalise-release/finalise-release.sh <release-issue-number>
```

For example, for "Finalise release 2246":

```bash
.claude/skills/finalise-release/finalise-release.sh 2246
```

The script prints progress to stderr and a summary to stdout: the Release
ticket URL and the issues moved (or the "nothing to move" notice). On any error
it exits non-zero with the reason — report that reason rather than retrying
blindly.

## Order of operations and safety

- The swimlane-existence check and the Release-ticket validation happen **before**
  any change, so a bad invocation aborts cleanly.
- All board moves happen **before** the Release ticket's description is edited.
  If a move fails (most likely a missing token scope), the run aborts *before*
  touching the ticket, leaving the board and the ticket consistent. Re-running
  after fixing the problem is safe: already-moved issues have left the
  `Approved...` swimlane and are simply skipped.
- The script never deletes anything and never touches issues outside the two
  swimlanes involved.

## Requirements

- `gh` CLI and `jq`.
- The active `gh` token needs:
  - **`project`** scope (write) to move items between swimlanes —
    `read:project` alone is **not** enough. If a move fails on permissions, refresh
    with `gh auth refresh -s project`.
  - rights to edit issues in `m2ms/fragalysis-frontend` (the repo is public, so
    `public_repo` suffices) to append to the Release ticket's description.

## Notes and assumptions

- "Swimlanes" are options of the board's single-select **Status** field; they are
  resolved by **name** at runtime, so changing their option IDs on the board does
  not break the skill (renaming them would — update the `*_SWIMLANE` constants in
  `finalise-release.sh` if they are ever renamed).
- `fragalysis-backend` is a real GitHub **issue label** (the same one used by the
  `create-backend-issue` skill), read from each issue's labels — not the board's
  "Target" field.
- The Release ticket and the moved issues both live in `m2ms/fragalysis-frontend`,
  so the appended bullets use same-repo `#<number>` references, which render as
  links with the title alongside (matching the existing Release-ticket style).
