import numpy as np
import pytest

from scripts.q2_relational_graph_qfi_gate import graph_fixtures, random_intertwiner
from scripts.q2_relational_gauge_quotient_gate import full_half_link_qfi
from scripts.q2_relational_joint_gauge_edge_quotient_gate import (
    analyze_sample,
    edge_singlet_common_generators,
    joint_quotient_basis,
    run,
)


def test_edge_common_generators_have_three_columns_per_internal_edge():
    for name, (nodes, edges) in graph_fixtures().items():
        common = edge_singlet_common_generators(nodes, edges)
        assert common.shape == (12 * nodes, 3 * len(edges))
        assert np.all(common.sum(axis=0) == 2.0)
        assert np.linalg.matrix_rank(common) == 3 * len(edges)


def test_joint_quotient_dimensions_for_path_and_cycle():
    for name, (nodes, edges) in graph_fixtures().items():
        quotient, info = joint_quotient_basis(nodes, edges)
        assert info["generator_families_independent"]
        assert info["joint_generator_rank"] == 3 * nodes + 3 * len(edges)
        assert quotient.shape[1] == 9 * nodes - 3 * len(edges)
        assert np.allclose(quotient.T @ quotient, np.eye(quotient.shape[1]), atol=1e-10)


def test_both_null_families_annihilate_full_half_link_qfi():
    rng = np.random.default_rng(20261118)
    for name, (nodes, edges) in graph_fixtures().items():
        states = [random_intertwiner(rng) for _ in range(nodes)]
        row = analyze_sample(states, edges)
        assert row["checks"]["vertex_gauss_generators_are_qfi_null"]
        assert row["checks"]["edge_singlet_common_generators_are_qfi_null"]
        assert row["checks"]["joint_quotient_has_no_remaining_null_modes"]


def test_joint_quotient_has_no_extra_projective_stabilizers():
    result = run(seed=20261118, trials=3)
    assert result["passed"]
    assert result["graphs"]["path_3"]["joint_quotient_dimension"] == 21
    assert result["graphs"]["cycle_4"]["joint_quotient_dimension"] == 24
    assert result["graphs"]["path_3"]["maximum_joint_quotient_nullity"] == 0
    assert result["graphs"]["cycle_4"]["maximum_joint_quotient_nullity"] == 0


def test_invalid_edge_common_map_is_rejected():
    with pytest.raises(ValueError):
        edge_singlet_common_generators(2, [(0, 0, 1, 1), (0, 0, 1, 2)])


def _closed_complete_five_graph(offset=0):
    edges = []
    for u in range(5):
        for v in range(u + 1, 5):
            lu = [x for x in range(5) if x != u].index(v)
            lv = [x for x in range(5) if x != v].index(u)
            edges.append((offset + u, lu, offset + v, lv))
    return edges


def test_topological_rank_theorem_for_open_and_closed_components():
    from scripts.q2_relational_joint_gauge_edge_quotient_gate import graph_topological_gauge_rank
    for nodes, edges in graph_fixtures().values():
        info = graph_topological_gauge_rank(nodes, edges)
        assert info["closed_components"] == 0
        assert info["topological_joint_generator_rank"] == 3 * (nodes + len(edges))
    info = graph_topological_gauge_rank(5, _closed_complete_five_graph())
    assert info["connected_components"] == 1
    assert info["closed_components"] == 1
    assert info["topological_generator_dependencies"] == 3
    assert info["topological_joint_generator_rank"] == 3 * (5 + 10 - 1)
    assert info["topological_joint_quotient_dimension"] == 18


def test_topological_rank_handles_disconnected_closed_components():
    from scripts.q2_relational_joint_gauge_edge_quotient_gate import graph_topological_gauge_rank
    edges = _closed_complete_five_graph(0) + _closed_complete_five_graph(5)
    info = graph_topological_gauge_rank(11, edges)
    assert info["connected_components"] == 3
    assert info["closed_components"] == 2
    assert info["open_components"] == 1
    assert info["topological_generator_dependencies"] == 6
    assert info["topological_joint_generator_rank"] == 3 * (11 + 20 - 2)


def test_joint_quotient_svd_rank_matches_topological_theorem():
    from scripts.q2_relational_joint_gauge_edge_quotient_gate import graph_topological_gauge_rank
    cases = list(graph_fixtures().values()) + [(5, _closed_complete_five_graph())]
    for nodes, edges in cases:
        _basis, info = joint_quotient_basis(nodes, edges)
        topology = graph_topological_gauge_rank(nodes, edges)
        assert info["svd_rank_matches_topology"]
        assert info["joint_generator_rank"] == topology["topological_joint_generator_rank"]
        assert info["joint_quotient_dimension"] == topology["topological_joint_quotient_dimension"]
