#!/usr/bin/env bash
#
# Copy a project-board planning issue into the fragalysis-backend repository.
#
# Given a project-board issue number, this:
#   1. reads the issue from the source repo (where the board's planning
#      issues live),
#   2. verifies it carries the "fragalysis-backend" label,
#   3. creates a copy (same title, same body) in the backend repo,
#   4. comments on the source issue with a link to the new backend issue.
#
# It is idempotent: a hidden marker is written into the comment so a second
# run for the same issue refuses to create a duplicate.
#
# Usage: create-backend-issue.sh <project-issue-number>

set -euo pipefail

# --- Configuration ----------------------------------------------------------
# The board (https://github.com/orgs/m2ms/projects/2) aggregates planning
# issues that live in SOURCE_REPO. The "fragalysis-backend" label only exists
# in that repo, so an issue carrying it must live there.
readonly SOURCE_REPO="m2ms/fragalysis-frontend"
readonly DEST_REPO="xchem/fragalysis-backend"
readonly REQUIRED_LABEL="fragalysis-backend"
# Hidden HTML marker used to detect an issue that has already been imported.
readonly MARKER="<!-- fragalysis-backend-import -->"

# --- Argument validation ----------------------------------------------------
if [ "$#" -ne 1 ]; then
  echo "Usage: $(basename "$0") <project-issue-number>" >&2
  exit 2
fi

issue_number="$1"
if ! [[ "${issue_number}" =~ ^[0-9]+$ ]]; then
  echo "Error: issue number must be a positive integer (got '${issue_number}')." >&2
  exit 2
fi

# --- Dependencies -----------------------------------------------------------
for cmd in gh jq; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Error: required command '${cmd}' is not installed." >&2
    exit 3
  fi
done

# --- Read the source issue --------------------------------------------------
echo "Reading ${SOURCE_REPO}#${issue_number} ..." >&2
if ! issue_json="$(gh issue view "${issue_number}" --repo "${SOURCE_REPO}" \
    --json number,title,body,labels,url,comments 2>/dev/null)"; then
  echo "Error: could not read issue ${SOURCE_REPO}#${issue_number}." \
       "Does it exist, and do you have access?" >&2
  exit 1
fi

title="$(jq -r '.title' <<<"${issue_json}")"
body="$(jq -r '.body' <<<"${issue_json}")"
source_url="$(jq -r '.url' <<<"${issue_json}")"
has_label="$(jq -r --arg l "${REQUIRED_LABEL}" \
  '[.labels[].name] | index($l) != null' <<<"${issue_json}")"

# --- Verify the required label ----------------------------------------------
if [ "${has_label}" != "true" ]; then
  echo "Error: ${SOURCE_REPO}#${issue_number} does not have the" \
       "'${REQUIRED_LABEL}' label; nothing to copy. Aborting." >&2
  exit 1
fi

# --- Refuse to create a duplicate -------------------------------------------
already_imported="$(jq -r --arg m "${MARKER}" \
  '[.comments[]? | select(.body | contains($m))] | length' <<<"${issue_json}")"
if [ "${already_imported}" -gt 0 ]; then
  echo "Error: ${SOURCE_REPO}#${issue_number} has already been imported" \
       "(import marker found in its comments). Aborting to avoid a duplicate." >&2
  exit 1
fi

# --- Build the new issue body -----------------------------------------------
# The body is a faithful copy of the source, with a footer that records the
# origin (and gives GitHub a back-link to the source issue).
body_file="$(mktemp)"
trap 'rm -f "${body_file}"' EXIT
{
  printf '%s\n' "${body}"
  printf '\n\n---\n'
  printf '_Imported from the project board — source issue: %s_\n' "${source_url}"
} >"${body_file}"

# --- Create the backend issue -----------------------------------------------
echo "Creating issue in ${DEST_REPO} ..." >&2
new_url="$(gh issue create --repo "${DEST_REPO}" \
  --title "${title}" --body-file "${body_file}")"

if [ -z "${new_url}" ]; then
  echo "Error: 'gh issue create' did not return a URL; the backend issue may" \
       "not have been created. Not commenting on the source issue." >&2
  exit 1
fi

# --- Comment on the source issue --------------------------------------------
comment_body="$(printf '%s\nBackend issue created from this ticket: %s' \
  "${MARKER}" "${new_url}")"
echo "Commenting on ${SOURCE_REPO}#${issue_number} ..." >&2
gh issue comment "${issue_number}" --repo "${SOURCE_REPO}" --body "${comment_body}"

# --- Report -----------------------------------------------------------------
echo "Done."
echo "  Source: ${source_url}"
echo "  Backend issue: ${new_url}"
