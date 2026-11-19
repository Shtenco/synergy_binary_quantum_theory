import pytest

import numpy as np

from scripts.q2_vertex_fusion_graph_gram import (
    _positive_spectrum,
    _relative_singular_rank,
    analyze_graph_block_fusion_gram,
)


@pytest.mark.parametrize(
    "graph_name,two_j,expected_rank,expected_boundary_dimension,expected_boundary_legs,expected_cycle_rank",
    [
        ("star_4", 1, 16, 42, 10, 0),
        ("cycle_4", 1, 12, 14, 8, 1),
        ("cycle_4", 2, 61, 91, 8, 1),
    ],
)
def test_tree_and_loop_blocks_have_distinct_fusion_gram_structure(
    graph_name,
    two_j,
    expected_rank,
    expected_boundary_dimension,
    expected_boundary_legs,
    expected_cycle_rank,
):
    result = analyze_graph_block_fusion_gram(graph_name, two_j)

    assert result["passed"], result
    assert result["graph"] == graph_name
    assert result["two_j"] == two_j
    assert result["cycle_rank"] == expected_cycle_rank
    assert result["boundary_leg_count"] == expected_boundary_legs
    assert result["boundary_invariant_dimension"] == expected_boundary_dimension
    assert result["fusion_map_rank"] == expected_rank
    assert result["fusion_map_rank_defect"] == result["input_dimension"] - expected_rank
    assert result["max_fused_boundary_gauss_relative_residual"] < 2e-9
    assert result["gram_eigenvalue_min"] > -2e-9

    if expected_cycle_rank == 0:
        dimension = two_j + 1
        edge_count = result["edge_count"]
        assert result["gram_is_scalar_identity"]
        assert result["tree_expected_gram_scale"] == pytest.approx(
            dimension ** (-2 * edge_count), abs=2e-12
        )
        assert result["tree_expected_gram_max_abs_error"] < 2e-12
        assert result["singular_value_condition_number_on_support"] == pytest.approx(
            1.0, abs=2e-10
        )
    else:
        assert result["tree_expected_gram_scale"] is None
        assert not result["gram_is_scalar_identity"]
        assert result["fusion_map_rank"] < result["input_dimension"]
        assert result["singular_value_condition_number_on_support"] > 2.0


def test_unknown_graph_fixture_is_rejected():
    with pytest.raises(ValueError, match="unknown graph fixture"):
        analyze_graph_block_fusion_gram("not_a_graph", 1)


@pytest.mark.parametrize("scale", [1.0, 1e-12, 1e-80])
def test_gram_rank_and_positive_spectrum_are_scale_invariant(scale):
    singular_values = scale * np.array([4.0, 2.0, 0.5, 1e-14])
    eigenvalues = singular_values**2

    assert _relative_singular_rank(singular_values, 1e-10) == 3
    positive = _positive_spectrum(eigenvalues, 1e-10)
    assert positive.size == 3
    assert np.all(positive > 0.0)
