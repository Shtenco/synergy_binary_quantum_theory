import numpy as np

from scripts.q2_relational_graph_qfi_gate import graph_fixtures, random_intertwiner
from scripts.q2_relational_two_sided_endpoint_qfi_gate import (
    analyze_sample,
    lifted_two_sided_qfi,
    two_sided_endpoint_injection,
    run,
)


def test_two_sided_injections_have_disjoint_endpoint_support():
    for name, (nodes, edges) in graph_fixtures().items():
        p_u, p_v, p_uv = two_sided_endpoint_injection(nodes, edges)
        assert p_u.shape == (12 * nodes, 3 * len(edges))
        assert p_v.shape == p_u.shape
        assert p_uv.shape == (12 * nodes, 6 * len(edges))
        assert np.all(p_u.sum(axis=0) == 1.0)
        assert np.all(p_v.sum(axis=0) == 1.0)
        assert np.all(p_u.T @ p_v == 0.0)


def test_two_sided_qfi_is_signed_lift_of_edge_qfi():
    rng = np.random.default_rng(20261118)
    for name, (nodes, edges) in graph_fixtures().items():
        row = analyze_sample([random_intertwiner(rng) for _ in range(nodes)], edges)
        assert all(row["checks"].values())


def test_common_endpoint_rotations_are_exact_null_directions():
    rng = np.random.default_rng(73)
    nodes, edges = graph_fixtures()["cycle_4"]
    row = analyze_sample([random_intertwiner(rng) for _ in range(nodes)], edges)
    assert row["two_sided_qfi_nullity"] >= 3 * len(edges)
    assert row["max_endpoint_tangent_sum_residual"] < 3e-9


def test_end_to_end_two_sided_endpoint_gate():
    result = run(seed=20261118, trials=2)
    assert result["passed"]
