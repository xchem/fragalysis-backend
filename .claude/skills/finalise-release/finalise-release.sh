#!/usr/bin/env bash
#
# Finalise a Release ticket on the m2ms Fragalysis project board.
#
# Given a "Release ..." ticket that lives in the board's "Releases" swimlane,
# this:
#   1. resolves the board's Status field and verifies the three required
#      swimlanes exist (Releases / Approved in staging - push to production /
#      In production (Done)) BEFORE doing any work,
#   2. validates the given issue is in "Releases" and its title starts with
#      "Release",
#   3. finds every issue in "Approved in staging - push to production" that
#      carries the "fragalysis-backend" label,
#   4. moves each of those issues to "In production (Done)",
#   5. appends a bullet-list (issue IDs + titles, as links) of the moved
#      issues to the Release ticket's description.
#
# If there is nothing to move it says so loudly and changes nothing — that may
# indicate a problem (e.g. the wrong Release ticket, or the issues were never
# approved), so the caller should relay it to the user.
#
# Usage: finalise-release.sh <release-issue-number>

set -euo pipefail

# --- Configuration ----------------------------------------------------------
# The board (https://github.com/orgs/m2ms/projects/2) is an org-level
# Projects-v2 board. Its planning issues — including the "Release ..." tickets —
# live in SOURCE_REPO. "Swimlanes" are options of the single-select "Status"
# field.
readonly ORG="m2ms"
readonly PROJECT_NUMBER="2"
readonly SOURCE_REPO="m2ms/fragalysis-frontend"
readonly STATUS_FIELD_NAME="Status"
readonly REQUIRED_LABEL="fragalysis-backend"

# The swimlanes this skill depends on. They are resolved by NAME at runtime so
# the script keeps working if the board's option IDs ever change.
readonly RELEASES_SWIMLANE="Releases"
readonly FROM_SWIMLANE="Approved in staging - push to production"
readonly TO_SWIMLANE="In production (Done)"

# --- Argument validation ----------------------------------------------------
if [ "$#" -ne 1 ]; then
  echo "Usage: $(basename "$0") <release-issue-number>" >&2
  exit 2
fi

release_number="$1"
if ! [[ "${release_number}" =~ ^[0-9]+$ ]]; then
  echo "Error: issue number must be a positive integer (got '${release_number}')." >&2
  exit 2
fi

# --- Dependencies -----------------------------------------------------------
for cmd in gh jq; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Error: required command '${cmd}' is not installed." >&2
    exit 3
  fi
done

# --- Resolve the board, the Status field and the required swimlanes ---------
# A single read query gives us the project node id, the Status field id, and
# every Status option (swimlane) with its id. We then check the three swimlanes
# we need all exist before touching anything.
echo "Resolving the project board and its Status field ..." >&2
fields_json="$(gh api graphql -f query='
  query($org: String!, $number: Int!) {
    organization(login: $org) {
      projectV2(number: $number) {
        id
        field(name: "'"${STATUS_FIELD_NAME}"'") {
          ... on ProjectV2SingleSelectField {
            id
            options { id name }
          }
        }
      }
    }
  }' -f org="${ORG}" -F number="${PROJECT_NUMBER}")"

project_id="$(jq -r '.data.organization.projectV2.id // empty' <<<"${fields_json}")"
status_field_id="$(jq -r '.data.organization.projectV2.field.id // empty' <<<"${fields_json}")"

if [ -z "${project_id}" ] || [ -z "${status_field_id}" ]; then
  echo "Error: could not read project ${ORG}/projects/${PROJECT_NUMBER} or its" \
       "'${STATUS_FIELD_NAME}' field. Check the board exists and your token has" \
       "the 'read:project' scope." >&2
  exit 1
fi

# Resolve each required swimlane to its option id, failing if any is missing.
option_id_for() {
  local wanted="$1"
  jq -r --arg n "${wanted}" \
    '.data.organization.projectV2.field.options[] | select(.name == $n) | .id' \
    <<<"${fields_json}"
}

releases_option_id="$(option_id_for "${RELEASES_SWIMLANE}")"
from_option_id="$(option_id_for "${FROM_SWIMLANE}")"
to_option_id="$(option_id_for "${TO_SWIMLANE}")"

missing=""
if [ -z "${releases_option_id}" ]; then missing="${missing}\n  - '${RELEASES_SWIMLANE}'"; fi
if [ -z "${from_option_id}" ];     then missing="${missing}\n  - '${FROM_SWIMLANE}'"; fi
if [ -z "${to_option_id}" ];       then missing="${missing}\n  - '${TO_SWIMLANE}'"; fi
if [ -n "${missing}" ]; then
  echo -e "Error: the board is missing required swimlane(s):${missing}" >&2
  echo "Aborting before any changes were made." >&2
  exit 1
fi
echo "All three required swimlanes exist." >&2

# --- Fetch every board item (paginated) -------------------------------------
# We need, for each item: its node id, its current Status, and its underlying
# issue's number/title/url/repo/labels. --paginate walks every page; jq -s
# slurps the per-page objects into one array we can query.
echo "Reading board items ..." >&2
items_json="$(gh api graphql --paginate -f query='
  query($org: String!, $number: Int!, $endCursor: String) {
    organization(login: $org) {
      projectV2(number: $number) {
        items(first: 100, after: $endCursor) {
          pageInfo { hasNextPage endCursor }
          nodes {
            id
            fieldValueByName(name: "'"${STATUS_FIELD_NAME}"'") {
              ... on ProjectV2ItemFieldSingleSelectValue { name }
            }
            content {
              ... on Issue {
                number
                title
                url
                repository { nameWithOwner }
                labels(first: 30) { nodes { name } }
              }
            }
          }
        }
      }
    }
  }' -f org="${ORG}" -F number="${PROJECT_NUMBER}" \
  | jq -s '[ .[].data.organization.projectV2.items.nodes[] ]')"

# --- Validate the Release ticket --------------------------------------------
# Locate the board item for the given issue number in the source repo.
release_item="$(jq -c --arg num "${release_number}" --arg repo "${SOURCE_REPO}" '
  map(select(.content.number == ($num | tonumber)
             and .content.repository.nameWithOwner == $repo))
  | first // empty' <<<"${items_json}")"

if [ -z "${release_item}" ]; then
  echo "Error: issue ${SOURCE_REPO}#${release_number} is not on the board." \
       "Pass the number of a 'Release ...' ticket from the '${RELEASES_SWIMLANE}'" \
       "swimlane." >&2
  exit 1
fi

release_status="$(jq -r '.fieldValueByName.name // ""' <<<"${release_item}")"
release_title="$(jq -r '.content.title // ""'        <<<"${release_item}")"
release_url="$(jq -r '.content.url // ""'            <<<"${release_item}")"

if [ "${release_status}" != "${RELEASES_SWIMLANE}" ]; then
  echo "Error: ${SOURCE_REPO}#${release_number} is in the '${release_status}'" \
       "swimlane, not '${RELEASES_SWIMLANE}'. Refusing to proceed." >&2
  exit 1
fi

# Title must start with the word "Release" (leading space tolerated, case
# insensitive).
release_title_lc="$(printf '%s' "${release_title}" | tr '[:upper:]' '[:lower:]' | sed -e 's/^[[:space:]]*//')"
if [[ "${release_title_lc}" != release* ]]; then
  echo "Error: ${SOURCE_REPO}#${release_number} title does not start with" \
       "'Release' (title: '${release_title}'). Refusing to proceed." >&2
  exit 1
fi
echo "Release ticket validated: #${release_number} '${release_title}'." >&2

# --- Find the issues to move ------------------------------------------------
# Items in the FROM swimlane whose underlying issue carries REQUIRED_LABEL.
# Emit one TSV line per candidate: <itemId>\t<number>\t<url>\t<title>.
candidates="$(jq -r \
  --arg from "${FROM_SWIMLANE}" \
  --arg label "${REQUIRED_LABEL}" '
  .[]
  | select(.fieldValueByName.name == $from)
  | select((.content.labels.nodes // []) | map(.name) | index($label) != null)
  | [.id, (.content.number|tostring), .content.url, .content.title]
  | @tsv' <<<"${items_json}")"

if [ -z "${candidates}" ]; then
  echo
  echo "NOTE: there are no issues in the '${FROM_SWIMLANE}' swimlane with the" \
       "'${REQUIRED_LABEL}' label, so nothing was moved and the Release ticket" \
       "was not changed."
  echo "This may indicate a problem — e.g. the issues were never approved into" \
       "'${FROM_SWIMLANE}', or you finalised the wrong Release ticket. Please" \
       "check the board: https://github.com/orgs/${ORG}/projects/${PROJECT_NUMBER}"
  exit 0
fi

num_candidates="$(printf '%s\n' "${candidates}" | wc -l | tr -d ' ')"
echo "Found ${num_candidates} issue(s) to move to '${TO_SWIMLANE}'." >&2

# --- Move each candidate to the TO swimlane ---------------------------------
# Mutating the board needs a token with the 'project' (write) scope. The first
# failed mutation aborts the run (set -e) before the Release body is touched, so
# the board and the ticket stay consistent and the run can be retried.
moved_bullets=""   # accumulates the markdown bullet list of moved issues
while IFS=$'\t' read -r item_id number url title; do
  [ -z "${item_id}" ] && continue
  echo "  moving #${number} (${title}) ..." >&2
  if ! mutation_out="$(gh api graphql -f query='
      mutation($project: ID!, $item: ID!, $field: ID!, $option: String!) {
        updateProjectV2ItemFieldValue(input: {
          projectId: $project
          itemId: $item
          fieldId: $field
          value: { singleSelectOptionId: $option }
        }) { projectV2Item { id } }
      }' \
      -f project="${project_id}" \
      -f item="${item_id}" \
      -f field="${status_field_id}" \
      -f option="${to_option_id}" 2>&1)"; then
    echo "Error: failed to move #${number} to '${TO_SWIMLANE}'." >&2
    echo "${mutation_out}" >&2
    echo >&2
    echo "If this is a permissions error, your token needs the 'project' (write)" \
         "scope. Refresh it with: gh auth refresh -s project" >&2
    echo "No change was made to the Release ticket. Re-run after fixing the" \
         "scope; already-moved issues will simply be skipped." >&2
    exit 1
  fi
  # Same-repo issue reference (#N) renders as a link with the title alongside.
  moved_bullets="${moved_bullets}- #${number} ${title}"$'\n'
done <<<"${candidates}"

# --- Append the moved issues to the Release ticket's description ------------
echo "Updating the Release ticket description ..." >&2
current_body="$(gh issue view "${release_number}" --repo "${SOURCE_REPO}" --json body -q '.body')"

body_file="$(mktemp)"
trap 'rm -f "${body_file}"' EXIT
{
  printf '%s\n' "${current_body}"
  printf '\n---\n'
  printf '### Backend issues released to production\n\n'
  printf '%s' "${moved_bullets}"
} >"${body_file}"

gh issue edit "${release_number}" --repo "${SOURCE_REPO}" --body-file "${body_file}" >/dev/null

# --- Report -----------------------------------------------------------------
echo
echo "Done. Moved ${num_candidates} issue(s) from '${FROM_SWIMLANE}' to" \
     "'${TO_SWIMLANE}' and listed them on the Release ticket:"
echo "  Release ticket: ${release_url}"
printf '%s' "${moved_bullets}" | sed 's/^- /    - /'
