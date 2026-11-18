import numpy as np

from scripts.q2_relational_graph_qfi_gate import graph_fixtures, random_intertwiner
from scripts.q2_relational_holonomy_boundary_tangent_gate import (
    compare_one,
    edge_probe_injection,
    run,
)


def test_edge_probe_injection_selects_first_endpoint_triplets():
    nodes, edges = graph_fixtures()["path_3"]
    mapping = edge_probe_injection(nodes, edges)
    assert mapping.shape == (12 * nodes, 3 * len(edges))
    assert np.all(mapping.sum(axis=0) == 1.0)
    assert np.linalg.matrix_rank(mapping) == 3 * len(edges)


def test_pullback_identity_on_path_and_cycle():
    rng = np.random.default_rng(20261118)
    for name in ("path_3", "cycle_4"):
        nodes, edges = graph_fixtures()[name]
        row = compare_one([random_intertwiner(rng) for _ in range(nodes)], edges)
        assert all(row["checks"].values())
        assert row["edge_projective_tangent_rank"] == row["edge_qfi_pullback_rank"]


def test_boundary_tangent_gate_end_to_end():
    result = run(seed=19, trials=2)
    assert result["passed"]
    assert all(
        graph["checks"]["direct_edge_qfi_matches_half_link_pullback"]
        for graph in result["graphs"].values()
    )
