import pytest

from scripts.q2_pl_generation3_stream import run


@pytest.fixture(scope="module")
def generation3_result():
    return run()


def test_generation3_direct_dual_and_causal_fibres_pass(generation3_result):
    result = generation3_result
    assert result["passed"] is True, result["checks"]
    assert all(result["checks"].values())


def test_generation3_matches_incidence_reference_at_smaller_levels(generation3_result):
    comparisons = generation3_result["small_level_incidence_comparisons"]
    assert [row["generation"] for row in comparisons] == [1, 2]
    assert all(row["exact_edge_and_route_match"] for row in comparisons)
    assert all(row["missing_from_direct"] == 0 for row in comparisons)
    assert all(row["extra_in_direct"] == 0 for row in comparisons)


def test_generation3_exact_sparse_graph_counts(generation3_result):
    counts = generation3_result["counts"]
    assert counts["parent_tetrahedra"] == 9_216
    assert counts["child_tetrahedra"] == 221_184
    assert counts["dual_edges"] == 442_368
    assert counts["unique_dual_edges"] == 442_368
    assert counts["route_edge_counts"] == [110_592] * 4
    assert counts["within_parent_edges"] == 331_776
    assert counts["across_parent_edges"] == 110_592
    assert counts["connected_components"] == 1


def test_generation3_causal_pushforward_is_surjective_and_conserved(generation3_result):
    counts = generation3_result["counts"]
    checks = generation3_result["checks"]
    assert counts["causal_fibres"] == 512
    assert counts["fibre_size_min"] > 0
    assert counts["fibre_size_max"] >= counts["fibre_size_min"]
    assert checks["all_causal_fibres_nonempty"] is True
    assert checks["pushforward_parent_mass_conserved"] is True
