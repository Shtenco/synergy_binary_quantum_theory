import pytest

pytest.importorskip("scipy")

from scripts.q2_causal_pl_correspondence_gate import run


@pytest.fixture(scope="module")
def feshbach_result():
    return run(max_generation=2, compute_feshbach=True)


def test_feshbach_correction_is_computed_and_positive_on_finite_refinements(feshbach_result):
    level0, level1, level2 = [row["feshbach"] for row in feshbach_result["levels"]]

    assert level0["leakage_frobenius_norm"] == 0.0
    assert level0["feshbach_correction_operator_norm"] == 0.0

    for level in (level1, level2):
        assert level["leakage_frobenius_norm"] > 0.0
        assert level["leakage_operator_norm"] > 0.0
        assert level["feshbach_correction_operator_norm"] > 0.0
        assert level["feshbach_correction_min_eigenvalue"] >= -1e-9
        assert level["effective_operator_min_eigenvalue"] >= -1e-8
        assert level["all_linear_solves_converged"] is True
        assert level["max_linear_solve_residual"] < 1e-8
        assert level["q_sector_gap_estimate"] is not None
        assert level["q_plus_p_largest_eigenvalue"] <= 8.0 + 1e-8


def test_generation_one_matches_exact_finite_control_values(feshbach_result):
    level1 = feshbach_result["levels"][1]["feshbach"]
    assert abs(level1["leakage_frobenius_norm"] ** 2 - 80.0 / 9.0) < 1e-9
    assert abs(level1["feshbach_correction_operator_norm"] - 40.0 / 51.0) < 1e-9
    assert abs(level1["q_sector_gap_estimate"] - 0.1522409349774257) < 1e-8
    assert abs(level1["feshbach_to_galerkin_operator_ratio"] - 20.0 / 119.0) < 1e-9


def test_generation_two_feshbach_correction_is_reproduced(feshbach_result):
    level2 = feshbach_result["levels"][2]["feshbach"]
    assert abs(level2["leakage_frobenius_norm"] - 9.98895305693888) < 1e-8
    assert abs(level2["feshbach_correction_operator_norm"] - 1.7366734310) < 1e-7
    assert abs(level2["q_sector_gap_estimate"] - 0.01171906316138909) < 1e-8
    assert abs(level2["feshbach_to_galerkin_operator_ratio"] - 0.3845362122) < 1e-7


def test_feshbach_correction_obeys_the_operator_lower_bound(feshbach_result):
    for row in feshbach_result["levels"][1:]:
        level = row["feshbach"]
        assert level["feshbach_correction_operator_norm"] >= (
            level["leakage_operator_norm"] ** 2 / 8.0 - 1e-8
        )
        assert level["feshbach_to_galerkin_operator_ratio"] > 0.0
