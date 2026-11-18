import numpy as np
import pytest

from scripts.q2_relational_graph_qfi_gate import (
    graph_fixtures,
    graph_qfi,
    random_intertwiner,
    run,
    validate_graph,
)


def test_path_and_cycle_full_graph_qfi_are_psd_with_cross_link_terms():
    result = run(seed=20261118, trials=8)
    assert result["passed"] is True
    for row in result["graphs"].values():
        assert row["minimum_qfi_eigenvalue"] >= -1e-10
        assert row["cross_link_qfi_frobenius_norm_max"] > 1e-3
    # The path reproduces the isolated two-link identity blocks; do not
    # assume this survives additional contractions around a cycle.
    assert result["graphs"]["path_3"]["diagonal_link_blocks_equal_identity"] is True
    # Closing the path into a cycle changes the boundary reduction: the
    # diagonal edge blocks are no longer universally the isolated I_3.
    assert result["graphs"]["cycle_4"]["max_diagonal_link_block_error_from_identity"] > 1e-5


def test_cycle_graph_has_eight_boundary_qubits_and_twelve_qfi_parameters():
    node_count, edges = graph_fixtures()["cycle_4"]
    assert node_count == 4
    assert len(edges) == 4
    rng = np.random.default_rng(7)
    states = [random_intertwiner(rng) for _ in range(node_count)]
    probability, qfi, symmetry, boundary = graph_qfi(states, edges)
    assert len(boundary) == 8
    assert qfi.shape == (12, 12)
    assert probability > 0
    assert symmetry < 1e-10


def test_graph_validation_rejects_reused_half_link():
    rng = np.random.default_rng(8)
    states = [random_intertwiner(rng) for _ in range(2)]
    with pytest.raises(ValueError, match="more than one graph edge"):
        validate_graph(states, [(0, 0, 1, 0), (0, 0, 1, 1)])


def test_graph_validation_rejects_non_normalized_node_state():
    states = [np.ones(16, dtype=complex), np.ones(16, dtype=complex)]
    with pytest.raises(ValueError, match="must be normalized"):
        validate_graph(states, [(0, 0, 1, 0)])
