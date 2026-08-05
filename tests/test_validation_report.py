#!/usr/bin/env python3
"""Tests for schema validation, metrics, classification and report outputs."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.generate_curvature_validation_report import (  # noqa: E402
    DataError,
    analyse,
    classify,
    classify_neutral,
    main,
    weighted_metrics,
)


def minimal_data() -> dict:
    return {
        "schema_version": 1,
        "metadata": {"title": "test", "recommendation": "not validated"},
        "implementation": {"summary": "test", "constants": [], "code_map": [], "risks": []},
        "benchmarks": [
            {
                "id": "case",
                "name": "case",
                "category": "neutral",
                "converged": True,
                "grid_independent": True,
                "datasets": [
                    {
                        "id": "q",
                        "quantity": "q",
                        "x": [0, 1],
                        "baseline": [2, 2],
                        "corrected": [1.1, 0.9],
                        "reference": [1, 1],
                        "weight": [1, 1],
                    }
                ],
            }
        ],
        "acceptance_criteria": [{"id": "unit", "criterion": "unit", "critical": True, "status": "pass"}],
    }


class MetricTests(unittest.TestCase):
    def test_weighted_metrics(self) -> None:
        metrics = weighted_metrics([2, 4], [1, 2], [1, 3])
        self.assertAlmostEqual(metrics.rmse, (13 / 4) ** 0.5)
        self.assertAlmostEqual(metrics.mae, 1.75)
        self.assertAlmostEqual(metrics.mean_bias, 1.75)
        self.assertEqual(metrics.maximum_absolute_error, 2)

    def test_classification(self) -> None:
        baseline = weighted_metrics([2], [1], [1])
        corrected = weighted_metrics([1.1], [1], [1])
        self.assertEqual(classify(baseline, corrected, converged=True, grid_independent=True, neutral_tolerance=0.02), "improves")
        self.assertEqual(classify(baseline, corrected, converged=False, grid_independent=True, neutral_tolerance=0.02), "non-conclusive")

    def test_neutral_classification_uses_baseline_change(self) -> None:
        self.assertEqual(
            classify_neutral(
                [1.0, 2.0],
                [1.01, 1.98],
                converged=True,
                grid_independent=None,
                neutral_tolerance=0.02,
            ),
            "neutral",
        )
        self.assertEqual(
            classify_neutral(
                [1.0],
                [1.03],
                converged=True,
                grid_independent=True,
                neutral_tolerance=0.02,
            ),
            "worsens",
        )

    def test_nonfinite_data_rejected(self) -> None:
        with self.assertRaises(DataError):
            weighted_metrics([float("nan")], [1], [1])


class GeneratorTests(unittest.TestCase):
    def test_generates_all_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "input.json"
            source.write_text(json.dumps(minimal_data()))
            output = directory / "report"
            self.assertEqual(main([str(source), "--output-dir", str(output)]), 0)
            self.assertTrue((output / "curvature-validation-report.md").is_file())
            self.assertTrue((output / "curvature-validation-report.html").is_file())
            self.assertTrue((output / "aggregated-metrics.json").is_file())
            self.assertTrue((output / "plots/case-q.svg").is_file())

    def test_csv_input_and_missing_reference_list(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "values.csv").write_text("x,baseline,corrected\n0,1,2\n")
            data = minimal_data()
            data["benchmarks"][0]["datasets"] = [{"id": "q", "csv": "values.csv"}]
            source = directory / "input.json"
            source.write_text(json.dumps(data))
            result = analyse(data, source, directory / "output")
            self.assertEqual(result["missing_reference"], ["case/q"])
            self.assertEqual(result["benchmarks"][0]["datasets"][0]["classification"], "non-conclusive")

    def test_schema_error_returns_two(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "bad.json"
            source.write_text("{}")
            self.assertEqual(main([str(source), "--output-dir", str(Path(temporary) / "out")]), 2)

    def test_critical_failure_returns_one_after_writing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            data = minimal_data()
            data["acceptance_criteria"][0]["status"] = "fail"
            source = directory / "input.json"
            source.write_text(json.dumps(data))
            output = directory / "out"
            self.assertEqual(main([str(source), "--output-dir", str(output)]), 1)
            self.assertTrue((output / "aggregated-metrics.json").is_file())

    def test_not_verified_only_blocks_strict_acceptance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            data = minimal_data()
            data["acceptance_criteria"][0]["status"] = "not-verified"
            source = directory / "input.json"
            source.write_text(json.dumps(data))
            normal_output = directory / "normal"
            strict_output = directory / "strict"
            self.assertEqual(
                main([str(source), "--output-dir", str(normal_output)]),
                0,
            )
            self.assertEqual(
                main(
                    [
                        str(source),
                        "--output-dir",
                        str(strict_output),
                        "--strict-acceptance",
                    ]
                ),
                1,
            )
            result = json.loads(
                (strict_output / "aggregated-metrics.json").read_text()
            )
            self.assertEqual(result["critical_incomplete"], ["unit"])


if __name__ == "__main__":
    unittest.main()
