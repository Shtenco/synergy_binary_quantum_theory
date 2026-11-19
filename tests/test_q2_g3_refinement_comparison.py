import json

import pytest

from scripts.q2_g3_refinement_comparison import compare_reports


def _write(path, data):
    path.write_text(json.dumps(data), encoding="utf-8")


def _fixture():
    times = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
    levels = []
    for generation in range(3):
        levels.append({
            "generation": generation,
            "operator_dimension": [1, 8, 64][generation],
            "operator_sparse_nnz": [0, 40, 424][generation],
            "operator_spectral_gap": [None, 1.0, 0.1][generation],
            "operator_minimum_eigenvalue": 0.0,
            "operator_largest_eigenvalue": 2.0 + generation,
            "operator_trace": 3.0 + generation,
            "operator_zero_mode_count": 1,
            "pushforward_measure_total": 16.0,
            "heat_trace": [
                {"time": t, "normalized_heat_trace": 1.0 / (1.0 + generation + t)}
                for t in times
            ],
            "effective_spectral_dimension": [
                {"time": t, "effective_spectral_dimension": 0.5 + generation * 0.1}
                for t in times[1:-1]
            ],
        })
    g3 = {
        "passed": True,
        "generation": 3,
        "sparse_galerkin": {
            "dimension": 512,
            "spectral_gap": 0.04,
            "zero_mode_count": 1,
            "measure_total": 16.0,
            "minimum_eigenvalue": 0.0,
            "maximum_eigenvalue": 3.0,
            "trace": 10.0,
            "nonzero_count": 4,
            "normalized_heat_trace": [
                {"t": t, "value": 1.0 / (4.0 + t)} for t in times
            ],
            "heat_trace_normalized": [
                {"t": t, "spectral_dimension": 0.8} for t in times[1:-1]
            ],
        },
    }
    return {"passed": True, "levels": levels}, g3


def test_g0_g3_comparison_reports_heat_trace_and_dimension_deltas(tmp_path):
    refinement, g3 = _fixture()
    refinement_path, g3_path = tmp_path / "refinement.json", tmp_path / "g3.json"
    _write(refinement_path, refinement)
    _write(g3_path, g3)

    report = compare_reports(refinement_path, g3_path)

    assert report["passed"] is True
    assert report["levels_compared"] == [0, 1, 2, 3]
    assert len(report["comparisons_to_g3"]) == 3
    assert len(report["level_summaries"]) == 4
    assert report["level_summaries"][2]["spectral_gap_ratio_to_g3"] == pytest.approx(2.5)
    assert report["level_summaries"][3]["spectral_gap_ratio_to_g3"] == 1.0
    assert report["comparisons_to_g3"][2]["max_absolute_heat_trace_delta"] > 0
    assert report["comparisons_to_g3"][2]["max_absolute_effective_dimension_delta"] == pytest.approx(0.1)
    assert "not established" in report["stability_conclusion"]


def test_g0_g3_comparison_rejects_missing_refinement_level(tmp_path):
    refinement, g3 = _fixture()
    refinement["levels"] = refinement["levels"][:2]
    refinement_path, g3_path = tmp_path / "refinement.json", tmp_path / "g3.json"
    _write(refinement_path, refinement)
    _write(g3_path, g3)

    with pytest.raises(ValueError, match="exactly the audited G0-G2"):
        compare_reports(refinement_path, g3_path)


def test_g0_g3_comparison_rejects_nonpassing_source_report(tmp_path):
    refinement, g3 = _fixture()
    g3["passed"] = False
    refinement_path, g3_path = tmp_path / "refinement.json", tmp_path / "g3.json"
    _write(refinement_path, refinement)
    _write(g3_path, g3)

    with pytest.raises(ValueError, match="must have passed"):
        compare_reports(refinement_path, g3_path)
