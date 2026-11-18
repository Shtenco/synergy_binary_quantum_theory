import numpy as np
import pytest

from scripts.q2_vertex_fusion_coarse_graining import (
    exhaustive_graph_partition_control,
    compare_staged_and_direct_fusion,
    deterministic_case,
    fuse_vertex_blocks,
)
from scripts.q2_relational_higher_spin_graph_qfi_gate import graph_fixtures, invariant_singlet_tensor, random_intertwiner
from scripts.q2_relational_higher_spin_network_qfi_gate import four_valent_intertwiner_basis


@pytest.mark.parametrize(
    "graph,blocks",
    [
        ("path_3", ((0, 1), (2,))),
        ("path_4", ((0, 1), (2, 3))),
        ("path_5", ((0, 1), (2, 3, 4))),
        ("star_4", ((0, 1), (2, 3))),
        ("cycle_4", ((0, 1), (2, 3))),
    ],
)
@pytest.mark.parametrize("two_j", [1, 2])
def test_exact_vertex_fusion_preserves_boundary_state_and_retained_edge_qfi(graph, blocks, two_j):
    result = deterministic_case(graph, two_j=two_j, blocks=blocks)
    assert result["passed"], result


def test_singleton_partition_reproduces_all_micro_edge_parameters():
    result = deterministic_case("path_4", two_j=1)
    assert result["passed"], result
    assert result["coarse_parameter_dimension"] == 9
    assert result["retained_edge_indices"] == [0, 1, 2]


@pytest.mark.parametrize("blocks", [
    ((0,), (1,)),
    ((0, 1),),
    ((0, 1), (1, 2)),
    ((0, 1), (2,), ()),
])
def test_invalid_vertex_partitions_are_rejected(blocks):
    nodes, edges = graph_fixtures()["path_3"]
    basis, _ = four_valent_intertwiner_basis(1)
    rng = np.random.default_rng(17)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    with pytest.raises(ValueError, match="partition|non-empty"):
        fuse_vertex_blocks(states, edges, invariant_singlet_tensor(1), blocks)


def test_vertex_fusion_rejects_reused_half_link():
    nodes, _ = graph_fixtures()["path_3"]
    basis, _ = four_valent_intertwiner_basis(1)
    rng = np.random.default_rng(18)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    with pytest.raises(ValueError, match="more than once"):
        fuse_vertex_blocks(states, [(0, 0, 1, 0), (0, 0, 2, 1)],
                           invariant_singlet_tensor(1), ((0, 1), (2,)))

@pytest.mark.parametrize("two_j", [1, 2])
def test_two_stage_vertex_fusion_equals_direct_fusion(two_j):
    nodes, edges = graph_fixtures()["path_5"]
    basis, _ = four_valent_intertwiner_basis(two_j)
    rng = np.random.default_rng(20261119 + two_j)
    states = [random_intertwiner(rng, basis) for _ in range(nodes)]
    singlet = invariant_singlet_tensor(two_j)
    result = compare_staged_and_direct_fusion(
        states,
        edges,
        singlet,
        first_blocks=((0, 1), (2,), (3, 4)),
        second_blocks=((0, 1), (2,)),
        direct_blocks=((0, 1, 2), (3, 4)),
    )
    assert result["passed"], result
    assert result["staged_retained_edge_indices"] == [2]
    assert result["direct_retained_edge_indices"] == [2]

@pytest.mark.parametrize("graph_name,expected_partitions", [
    ("path_3", 5),   # Bell number B_3
    ("path_4", 15),  # Bell number B_4
    ("path_5", 52),  # Bell number B_5
    ("star_4", 15),
    ("cycle_4", 15),
])
@pytest.mark.parametrize("two_j", [1, 2])
def test_exhaustive_graph_vertex_partitions_and_block_order_invariance(
    graph_name, expected_partitions, two_j
):
    result = exhaustive_graph_partition_control(
        graph_name, two_j=two_j, seed=20261120 + two_j
    )
    assert result["graph"] == graph_name
    assert result["partition_count"] == expected_partitions
    assert result["expected_partition_count"] == expected_partitions
    assert result["includes_block_order_reversal"]
    assert result["max_fused_block_gauss_relative_residual"] < 2e-9
    assert result["passed"], result


