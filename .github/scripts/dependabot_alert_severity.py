#!/usr/bin/env python3
"""Decide if a Dependabot pull request on the 'security' branch must auto-merge.

Issue #1088. Dependabot sends *version*-update pull requests to the 'security'
branch. GitHub does not associate a version update with a security alert, so
'dependabot/fetch-metadata' reports a CVSS of 0 for them. This script instead
finds the open Dependabot alerts for the packages that the pull request updates,
and compares the highest alert severity with a minimum severity.

Inputs:
    --alerts        A JSON file that holds a list of Dependabot alerts, as the
                    GitHub REST API returns them
                    (GET /repos/{owner}/{repo}/dependabot/alerts).
    --dependencies  A comma-separated list of the updated package names
                    (the 'dependency-names' output of fetch-metadata).
    --min-severity  One of: low, medium, high, critical.

Output:
    'true' (merge) or 'false' (do not merge) on stdout.
    An invalid minimum severity is an error (non-zero exit).

The script uses only the Python standard library.
"""

import argparse
import json
import re
import sys

# GitHub advisory severities, lowest first.
SEVERITIES = ("low", "medium", "high", "critical")


def normalise_name(name):
    """Return the PEP 503 normalised form of a package name."""
    return re.sub(r"[-_.]+", "-", name.strip()).lower()


def highest_severity(alerts, dependency_names):
    """Return the highest severity of the open alerts for the named packages.

    Returns None if no open alert (with a known severity) matches.
    """
    names = {normalise_name(name) for name in dependency_names if name.strip()}
    highest = None
    for alert in alerts:
        if alert.get("state") != "open":
            continue
        package_name = alert.get("dependency", {}).get("package", {}).get("name", "")
        if normalise_name(package_name) not in names:
            continue
        severity = alert.get("security_advisory", {}).get("severity", "").lower()
        if severity not in SEVERITIES:
            continue
        if highest is None or SEVERITIES.index(severity) > SEVERITIES.index(highest):
            highest = severity
    return highest


def should_merge(alerts, dependency_names, min_severity):
    """Return True if an open alert is at or above the minimum severity.

    Raises ValueError if min_severity is not a known severity.
    """
    threshold = min_severity.strip().lower()
    if threshold not in SEVERITIES:
        raise ValueError(
            f"Invalid minimum severity '{min_severity}'."
            f" Use one of: {', '.join(SEVERITIES)}"
        )
    highest = highest_severity(alerts, dependency_names)
    if highest is None:
        return False
    return SEVERITIES.index(highest) >= SEVERITIES.index(threshold)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--alerts", required=True)
    parser.add_argument("--dependencies", required=True)
    parser.add_argument("--min-severity", required=True)
    args = parser.parse_args()

    with open(args.alerts, encoding="utf-8") as alerts_file:
        alerts = json.load(alerts_file)
    try:
        merge = should_merge(alerts, args.dependencies.split(","), args.min_severity)
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1
    print("true" if merge else "false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
