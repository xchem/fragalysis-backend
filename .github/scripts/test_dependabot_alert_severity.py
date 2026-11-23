"""Unit tests for dependabot_alert_severity.py (issue #1088).

Run with: python -m unittest discover -s .github/scripts
"""

import json
import os
import subprocess
import sys
import tempfile
import unittest

from dependabot_alert_severity import highest_severity, normalise_name, should_merge

_SCRIPT = os.path.join(os.path.dirname(__file__), "dependabot_alert_severity.py")


def _alert(name, severity, state="open", ecosystem="pip"):
    """Return a minimal Dependabot alert, as the GitHub alerts API returns it."""
    return {
        "state": state,
        "dependency": {"package": {"ecosystem": ecosystem, "name": name}},
        "security_advisory": {"severity": severity},
    }


class NormaliseNameTest(unittest.TestCase):
    def test_pep503_normalisation(self):
        self.assertEqual(
            normalise_name("Django_REST.framework"), "django-rest-framework"
        )

    def test_strips_whitespace(self):
        self.assertEqual(normalise_name(" requests "), "requests")


class HighestSeverityTest(unittest.TestCase):
    def test_no_alerts(self):
        self.assertIsNone(highest_severity([], ["django"]))

    def test_alert_for_another_package(self):
        self.assertIsNone(highest_severity([_alert("urllib3", "critical")], ["django"]))

    def test_matches_a_normalised_name(self):
        alerts = [_alert("Django", "high")]
        self.assertEqual(highest_severity(alerts, ["django"]), "high")

    def test_returns_the_highest_of_several(self):
        alerts = [
            _alert("django", "low"),
            _alert("django", "critical"),
            _alert("django", "medium"),
        ]
        self.assertEqual(highest_severity(alerts, ["django"]), "critical")

    def test_considers_every_dependency_of_a_grouped_pull_request(self):
        alerts = [_alert("django", "medium"), _alert("urllib3", "high")]
        self.assertEqual(highest_severity(alerts, ["django", "urllib3"]), "high")

    def test_ignores_alerts_that_are_not_open(self):
        alerts = [
            _alert("django", "critical", state="dismissed"),
            _alert("django", "critical", state="fixed"),
        ]
        self.assertIsNone(highest_severity(alerts, ["django"]))

    def test_ignores_an_unknown_severity(self):
        self.assertIsNone(highest_severity([_alert("django", "unknown")], ["django"]))


class ShouldMergeTest(unittest.TestCase):
    def test_severity_at_threshold_merges(self):
        self.assertTrue(should_merge([_alert("django", "high")], ["django"], "high"))

    def test_severity_above_threshold_merges(self):
        self.assertTrue(
            should_merge([_alert("django", "critical")], ["django"], "high")
        )

    def test_severity_below_threshold_does_not_merge(self):
        self.assertFalse(should_merge([_alert("django", "medium")], ["django"], "high"))

    def test_no_alert_does_not_merge(self):
        self.assertFalse(should_merge([], ["django"], "low"))

    def test_threshold_is_case_insensitive(self):
        self.assertTrue(should_merge([_alert("django", "high")], ["django"], "HIGH"))

    def test_invalid_threshold_is_an_error(self):
        with self.assertRaises(ValueError):
            should_merge([_alert("django", "high")], ["django"], "severe")


class CommandLineTest(unittest.TestCase):
    def _run(self, alerts, dependencies, min_severity):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(alerts, f)
        try:
            return subprocess.run(
                [
                    sys.executable,
                    _SCRIPT,
                    "--alerts",
                    f.name,
                    "--dependencies",
                    dependencies,
                    "--min-severity",
                    min_severity,
                ],
                capture_output=True,
                text=True,
                check=False,
            )
        finally:
            os.unlink(f.name)

    def test_prints_true_when_the_pull_request_must_merge(self):
        result = self._run([_alert("django", "high")], "django, urllib3", "high")
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "true")

    def test_prints_false_when_the_pull_request_must_not_merge(self):
        result = self._run([_alert("django", "low")], "django", "high")
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "false")

    def test_invalid_threshold_fails(self):
        result = self._run([], "django", "severe")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("severe", result.stderr)


if __name__ == "__main__":
    unittest.main()
