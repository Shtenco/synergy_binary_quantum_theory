import numpy as np
import pytest

from scripts.q2_vertex_fusion_haar_pushforward import analyze_independent_haar_fusion


@pytest.mark.parametrize(
    "two_j,expected_invariant_dimension",
    [(1, 5), (2, 15), (3, 34), (4, 65)],
)
def test_local_haar_pushforward_has_proper_six_valent_invariant_support(
    two_j, expected_invariant_dimension
):
    result = analyze_independent_haar_fusion(two_j)
    local_dimension = two_j + 1
    input_dimension = local_dimension**2
    assert result["local_intertwiner_dimension"] == local_dimension
    assert result["product_input_dimension"] == input_dimension
    assert result["fused_boundary_leg_count"] == 6
    assert result["full_six_valent_invariant_dimension"] == expected_invariant_dimension
    assert result["fusion_map_rank"] == input_dimension
    assert np.allclose(
        result["fusion_map_singular_values"],
        [1.0 / local_dimension] * input_dimension,
        atol=2e-9,
    )
    assert result["fusion_map_gram_max_abs_error_from_I_over_d_squared"] < 2e-9
    assert result["six_valent_gauss_max_relative_residual"] < 2e-9
    assert result["pushforward_covariance_rank"] == input_dimension
    assert np.allclose(
        result["pushforward_nonzero_covariance_eigenvalues"],
        [1.0 / input_dimension] * input_dimension,
        atol=2e-9,
    )
    assert result["pushforward_covariance_purity"] == pytest.approx(
        1.0 / input_dimension, abs=2e-9
    )
    assert result["full_six_valent_haar_covariance_purity"] == pytest.approx(
        1.0 / expected_invariant_dimension, abs=2e-9
    )
    expected_distance = np.sqrt(
        1.0 / input_dimension - 1.0 / expected_invariant_dimension
    )
    assert result["pushforward_vs_full_haar_covariance_frobenius_distance"] == pytest.approx(
        expected_distance, abs=2e-9
    )
    assert result["passed"], result
