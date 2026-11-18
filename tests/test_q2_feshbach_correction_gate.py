import pytest

pytest.importorskip("scipy")

from scripts.q2_causal_pl_correspondence_gate import run


def test_feshbach_correction_is_computed_and_positive_on_finite_refinements():
    result = run(max_generation=1, compute_feshbach=True)
    level0, level1 = result["levels"]
    f0 = level0["feshbach"]
    f1 = level1["feshbach"]

    assert f0["leakage_frobenius_norm"] == 0.0
    assert f0["feshbach_correction_operator_norm"] == 0.0
    assert f1["leakage_frobenius_norm"] > 0.0
    assert f1["leakage_operator_norm"] > 0.0
    assert f1["feshbach_correction_operator_norm"] > 0.0
    assert f1["feshbach_correction_min_eigenvalue"] >= -1e-9
    assert f1["effective_operator_min_eigenvalue"] >= -1e-8
    assert f1["all_linear_solves_converged"] is True
    assert f1["max_linear_solve_residual"] < 1e-8
    assert abs(f1["leakage_frobenius_norm"] ** 2 - 80.0 / 9.0) < 1e-9
    assert abs(f1["feshbach_correction_operator_norm"] - 40.0 / 51.0) < 1e-9
    assert abs(f1["q_sector_gap_estimate"] - 0.1522409349774257) < 1e-8
    assert f1["q_plus_p_condition_estimate"] > 1.0


def test_feshbach_correction_obeys_the_operator_lower_bound():
    result = run(max_generation=1, compute_feshbach=True)
    f1 = result["levels"][1]["feshbach"]
    assert f1["feshbach_correction_operator_norm"] >= (
        f1["leakage_operator_norm"] ** 2 / 8.0 - 1e-9
    )
    assert f1["feshbach_to_galerkin_operator_ratio"] > 0.0
