import numpy as np
import pytest

from scripts.q2_vertex_fusion_haar_pushforward import analyze_independent_haar_fusion


def test_two_vertex_local_haar_pushforward_is_rank_four_inside_six_valent_invariants():
    result = analyze_independent_haar_fusion()
    assert result["local_intertwiner_dimension"] == 2
    assert result["product_input_dimension"] == 4
    assert result["fused_boundary_leg_count"] == 6
    assert result["full_six_valent_invariant_dimension"] == 5
    assert result["fusion_map_rank"] == 4
    assert np.allclose(result["fusion_map_singular_values"], [0.5] * 4, atol=2e-10)
    assert result["fusion_map_gram_max_abs_error_from_I_over_4"] < 2e-10
    assert result["fusion_image_invariant_subspace_residual"] < 2e-10
    assert result["pushforward_covariance_rank"] == 4
    assert np.allclose(
        result["pushforward_nonzero_covariance_eigenvalues"],
        [0.25] * 4,
        atol=2e-10,
    )
    assert result["pushforward_covariance_purity"] == pytest.approx(0.25, abs=2e-10)
    assert result["full_six_valent_haar_covariance_purity"] == pytest.approx(0.2, abs=2e-10)
    assert result["pushforward_vs_full_haar_covariance_frobenius_distance"] == pytest.approx(
        np.sqrt(0.05), abs=2e-10
    )
    assert result["passed"], result


def test_pushforward_report_does_not_claim_physical_bqg_measure():
    result = analyze_independent_haar_fusion()
    assert "not a statement about the unspecified physical BQG measure" in result["interpretation"]
    assert "only" in result["claim_boundary"]
