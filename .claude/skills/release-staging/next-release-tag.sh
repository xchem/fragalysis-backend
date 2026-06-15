#!/usr/bin/env bash
# Compute the next production release tag and the prior full-release tag.
#
# Prints exactly two lines to stdout, e.g.:
#   NEXT_TAG=2026.06.2
#   PRIOR_TAG=2026.06.1
#
# Tag scheme is YYYY.MM.ITERATION where ITERATION restarts at 1 each calendar
# month (UTC). NEXT_TAG is one greater than the highest full-release iteration
# already used this month, or .1 if this is the month's first release.
# PRIOR_TAG is the most recent full (non-prerelease) GitHub release, used as the
# starting point for auto-generated release notes. It is empty if none exist.
#
# Pre-release tags (e.g. "2026.06.1-rc.1") are ignored for both calculations.

set -euo pipefail

REPO="xchem/fragalysis-backend"

# Year.month in UTC, e.g. "2026.06". Dots are escaped for the regex below.
ym=$(date -u +%Y.%m)
ym_re="^${ym//./\\.}\.([0-9]+)$"

prior_tag=""
max_iter=0

# gh lists releases newest-first, so the first non-prerelease is the prior one.
while IFS= read -r tag; do
  if [ -z "${tag}" ]; then
    continue
  fi
  if [ -z "${prior_tag}" ]; then
    prior_tag="${tag}"
  fi
  if [[ "${tag}" =~ ${ym_re} ]]; then
    iter="${BASH_REMATCH[1]}"
    if [ "${iter}" -gt "${max_iter}" ]; then
      max_iter="${iter}"
    fi
  fi
done < <(gh release list --repo "${REPO}" --limit 200 \
           --json tagName,isPrerelease \
           --jq '.[] | select(.isPrerelease == false) | .tagName')

next_iter=$((max_iter + 1))

echo "NEXT_TAG=${ym}.${next_iter}"
echo "PRIOR_TAG=${prior_tag}"
