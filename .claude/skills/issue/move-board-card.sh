#!/usr/bin/env bash
#
# Move the project-board card for a backend issue to a target swimlane.
#
# The m2ms project board (https://github.com/orgs/m2ms/projects/2) tracks the
# *frontend* planning issues (m2ms/fragalysis-frontend), not the backend issues
# this skill works on. A backend issue created by `create-backend-issue` carries
# a footer linking to its source frontend issue; we follow that link to find the
# board card.
#
# Target (second argument) is one of:
#   --in-progress   Move to "{Category} - In Progress", where Category is the
#                   text before the first " - " in the card's current lane
#                   (e.g. "Infra - Backlog" -> "Infra - In Progress"). Cards that
#                   are *already* in an "In Progress" lane are left untouched.
#   <lane name>     Move to that exact lane, matched case-insensitively
#                   (e.g. "Dev Done - Do review (DEV)").
#
# If the target lane does not exist, no board card can be found, or the card is
# already in the target lane, nothing is moved and the reason is printed; none of
# those are fatal to the issue workflow, so they exit 0.
#
# Usage: move-board-card.sh <backend-issue-number> (--in-progress | <lane name>)

set -euo pipefail

# --- Configuration ----------------------------------------------------------
readonly PROJECT_NUMBER="2"
readonly FRONTEND_OWNER="m2ms"
readonly FRONTEND_REPO="fragalysis-frontend"
readonly BACKEND_REPO="xchem/fragalysis-backend"
readonly STATUS_FIELD_NAME="Status"
# Suffix appended to a category for the --in-progress target, matched
# case-insensitively (the board mixes "In Progress" and "In progress").
readonly IN_PROGRESS_SUFFIX="In Progress"

# Lowercase helper (the host's bash is 3.2, which lacks ${var,,}).
lc() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }

# --- Argument validation ----------------------------------------------------
if [ "$#" -ne 2 ]; then
  echo "Usage: $(basename "$0") <backend-issue-number> (--in-progress | <lane name>)" >&2
  exit 2
fi

backend_number="$1"
target_arg="$2"
if ! [[ "${backend_number}" =~ ^[0-9]+$ ]]; then
  echo "Error: issue number must be a positive integer (got '${backend_number}')." >&2
  exit 2
fi
if [ -z "${target_arg}" ]; then
  echo "Error: target lane (or --in-progress) must not be empty." >&2
  exit 2
fi

# --- Dependencies -----------------------------------------------------------
for cmd in gh jq; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Error: required command '${cmd}' is not installed." >&2
    exit 3
  fi
done

# --- Find the source (board) issue ------------------------------------------
# The board card lives on the frontend issue linked from the backend issue body.
if ! body="$(gh issue view "${backend_number}" --repo "${BACKEND_REPO}" \
    --json body -q '.body' 2>/dev/null)"; then
  echo "Error: could not read ${BACKEND_REPO}#${backend_number}." >&2
  exit 1
fi

frontend_number="$(grep -oiE "${FRONTEND_REPO}/issues/[0-9]+" <<<"${body}" \
  | head -1 | grep -oE '[0-9]+$' || true)"

if [ -z "${frontend_number}" ]; then
  echo "Couldn't move the board card: ${BACKEND_REPO}#${backend_number} has no" \
       "link to a ${FRONTEND_OWNER}/${FRONTEND_REPO} source issue, so there is no" \
       "board card to find."
  exit 0
fi

# --- Read the card's current lane and the board's Status options ------------
item_json="$(gh api graphql -f query='
  query($owner: String!, $repo: String!, $number: Int!) {
    repository(owner: $owner, name: $repo) {
      issue(number: $number) {
        projectItems(first: 20) {
          nodes {
            id
            project {
              number
              id
              field(name: "'"${STATUS_FIELD_NAME}"'") {
                ... on ProjectV2SingleSelectField {
                  id
                  options { id name }
                }
              }
            }
            fieldValueByName(name: "'"${STATUS_FIELD_NAME}"'") {
              ... on ProjectV2ItemFieldSingleSelectValue { name }
            }
          }
        }
      }
    }
  }' -f owner="${FRONTEND_OWNER}" -f repo="${FRONTEND_REPO}" \
     -F number="${frontend_number}")"

# The item on *our* board (project number PROJECT_NUMBER).
node="$(jq -c --argjson pn "${PROJECT_NUMBER}" \
  '.data.repository.issue.projectItems.nodes[]? | select(.project.number == $pn)' \
  <<<"${item_json}" | head -1)"

if [ -z "${node}" ]; then
  echo "Couldn't move the board card: ${FRONTEND_OWNER}/${FRONTEND_REPO}" \
       "#${frontend_number} is not on project board #${PROJECT_NUMBER}."
  exit 0
fi

item_id="$(jq -r '.id' <<<"${node}")"
project_id="$(jq -r '.project.id' <<<"${node}")"
status_field_id="$(jq -r '.project.field.id' <<<"${node}")"
current_lane="$(jq -r '.fieldValueByName.name // ""' <<<"${node}")"

if [ -z "${current_lane}" ]; then
  echo "Couldn't move the board card: ${FRONTEND_OWNER}/${FRONTEND_REPO}" \
       "#${frontend_number} has no Status set on the board."
  exit 0
fi

# --- Work out the target lane -----------------------------------------------
if [ "${target_arg}" = "--in-progress" ]; then
  # Only move cards that are not already in progress. In-progress lanes are named
  # "{Category} - In Progress" (or lowercase "In progress"), so treat any lane
  # whose name ends in "in progress" as already in progress.
  if [[ "$(lc "${current_lane}")" == *"in progress" ]]; then
    echo "Board card for #${frontend_number} is already in an 'In Progress' lane" \
         "('${current_lane}'). Not moving it."
    exit 0
  fi
  # Category = the text before the first " - " (the whole name if there is none).
  category="${current_lane%% - *}"
  target_name="${category} - ${IN_PROGRESS_SUFFIX}"
else
  target_name="${target_arg}"
fi

# Already in the target lane? Nothing to do.
if [ "$(lc "${current_lane}")" = "$(lc "${target_name}")" ]; then
  echo "Board card for #${frontend_number} is already in '${current_lane}'." \
       "Nothing to move."
  exit 0
fi

# Find the matching option id, comparing names case-insensitively.
target_option_id="$(jq -r --arg want "$(lc "${target_name}")" \
  '.project.field.options[] | select((.name | ascii_downcase) == $want) | .id' \
  <<<"${node}" | head -1)"

if [ -z "${target_option_id}" ]; then
  echo "Couldn't move the board card: no '${target_name}' lane exists on board" \
       "#${PROJECT_NUMBER} (current lane '${current_lane}'). Leaving it where it is."
  exit 0
fi

# --- Move it ----------------------------------------------------------------
echo "Moving board card for #${frontend_number} from '${current_lane}' to" \
     "'${target_name}' ..." >&2
if ! gh api graphql -f query='
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
    -f option="${target_option_id}" >/dev/null; then
  echo "Error: failed to move the board card. If this is a permissions error," \
       "your token needs the 'project' (write) scope (gh auth refresh -s project)." >&2
  exit 1
fi

echo "Done. Board card for ${FRONTEND_OWNER}/${FRONTEND_REPO}#${frontend_number}" \
     "moved to '${target_name}'."
