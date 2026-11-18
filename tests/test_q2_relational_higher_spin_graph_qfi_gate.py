import numpy as np
import pytest

from scripts.q2_relational_graph_qfi_gate import graph_fixtures, random_intertwiner as spin_half_intertwiner, graph_qfi as spin_half_graph_qfi
from scripts.q2_relational_higher_spin_network_qfi_gate import (
    four_valent_intertwiner_basis,
    invariant_singlet_tensor,
    random_intertwiner,
    spin_generators,
)
from scripts.q2_relational_higher_spin_graph_qfi_gate import (
    analyze_graph,
    contract_graph,
    graph_qfi,
    run,
)


@pytest.mark.parametrize("two_j", [1, 2])
def test_higher_spin_path_and_cycle_qfi(two_j):
    result = run(two_j_values=(two_j,), trials=2)
    assert result["passed"]


def test_spin_half_graph_generalization_matches_existing_graph_gate():
    rng = np.random.default_rng(20261118)
    basis, _ = four_valent_intertwiner_basis(1)
    singlet = invariant_singlet_tensor(1)
    generators = spin_generators(1)
    for name, (nodes, edges) in graph_fixtures().items():
        states = [random_intertwiner(rng, basis) for _ in range(nodes)]
        _p, qfi_new, _s, _b = graph_qfi(states, edges, singlet, generators)
        _p2, qfi_old, _s2, _b2 = spin_half_graph_qfi(states, edges)
        assert np.allclose(qfi_new, qfi_old, atol=2e-9, rtol=2e-9)


def test_invalid_graph_reuses_half_link_rejected():
    basis, _ = four_valent_intertwiner_basis(1)
    states = [random_intertwiner(np.random.default_rng(4), basis) for _ in range(2)]
    singlet = invariant_singlet_tensor(1)
    with pytest.raises(ValueError):
        contract_graph(states, [(0, 0, 1, 1), (0, 0, 1, 2)], singlet)

def test_extended_graph_family_reports_complete_qfi_spectral_partition():
    result = run(two_j_values=(1, 2), trials=2)
    assert result["passed"]
    assert {"path_3", "path_4", "path_5", "star_4", "cycle_4"} <= set(result["graphs"]["0.5"])
    expected_boundary_legs = {
        "path_3": 8,
        "path_4": 10,
        "path_5": 12,
        "star_4": 10,
        "cycle_4": 8,
    }
    for spin, graph_rows in result["graphs"].items():
        for name, expected_boundary in expected_boundary_legs.items():
            if name not in graph_rows:
                continue
            row = graph_rows[name]
            assert row["boundary_half_links"] == expected_boundary
            assert row["qfi_parameter_dimension"] == 3 * row["edges"]
            assert row["qfi_rank_max_over_trials"] <= row["qfi_parameter_dimension"]
            assert row["qfi_rank_min_over_trials"] >= 0
            assert row["qfi_rank_min_over_trials"] + row["qfi_nullity_min_over_trials"] <= row["qfi_parameter_dimension"]
            assert row["checks"]["qfi_rank_nullity_partition_is_consistent"]
            assert row["qfi_spectral_maximum_over_trials"] >= row["qfi_spectral_minimum_over_trials"]

@pytest.mark.parametrize(
    "edge",
    [
        (0.5, 0, 1, 1),
        (0, 0, 1, 1.0),
        (False, 0, 1, 1),
        (0, 0, 1, "1"),
    ],
)
def test_higher_spin_graph_rejects_non_integer_half_link_indices(edge):
    basis, _ = four_valent_intertwiner_basis(1)
    states = [random_intertwiner(np.random.default_rng(8), basis) for _ in range(2)]
    singlet = invariant_singlet_tensor(1)
    with pytest.raises(ValueError, match="indices must be integers"):
        contract_graph(states, [edge], singlet)



def test_path5_dimensionless_spectral_diagnostics_are_reported():
    result = run(two_j_values=(1,), trials=2)
    row = result["graphs"]["0.5"]["path_5"]
    assert row["nodes"] == 5
    assert row["edges"] == 4
    assert row["boundary_half_links"] == 12
    assert row["qfi_parameter_dimension"] == 12
    assert row["qfi_spectral_normalization"].startswith("eigenvalues divided")
    assert set(row["normalized_spectral_moments_mean"]) == {"1", "2", "3", "4"}
    assert len(row["normalized_spectral_quantiles_mean"]) == 5
    assert 0.0 <= row["qfi_rank_fraction_mean"] <= 1.0
    assert row["normalized_spectral_variance_mean"] >= -1e-10
    assert result["graph_spin_coverage"]["path_5"] == [0.5]


@pytest.mark.parametrize("two_j", [3, 4])
def test_high_spin_graph_spectra_on_bounded_boundary_fixtures(two_j):
    # Keep only 8-boundary-leg fixtures so j=3/2,2 remain tractable.
    result = run(
        two_j_values=(two_j,),
        trials=2,
        graph_names=("path_3", "cycle_4"),
    )
    assert result["passed"]
    spin_key = str(two_j / 2.0)
    assert set(result["graphs"][spin_key]) == {"path_3", "cycle_4"}
    for row in result["graphs"][spin_key].values():
        assert row["boundary_half_links"] == 8
        assert row["qfi_rank_fraction_mean"] <= 1.0
        assert row["checks"]["qfi_is_psd"]


def test_collective_edge_parameter_map_and_pullback_are_exact():
    from scripts.q2_relational_higher_spin_graph_qfi_gate import (
        collective_edge_parameter_map,
        collective_edge_pullback_qfi,
    )
    blocks = ((0, 1), (2,))
    mapping = collective_edge_parameter_map(3, blocks)
    expected = np.zeros((9, 6))
    expected[0:3, 0:3] = np.eye(3)
    expected[3:6, 0:3] = np.eye(3)
    expected[6:9, 3:6] = np.eye(3)
    assert np.array_equal(mapping, expected)
    rng = np.random.default_rng(20261118)
    raw = rng.normal(size=(9, 9))
    micro_qfi = raw.T @ raw
    effective, returned_map = collective_edge_pullback_qfi(micro_qfi, blocks)
    assert np.array_equal(returned_map, mapping)
    assert np.allclose(effective, mapping.T @ micro_qfi @ mapping)
    assert np.linalg.eigvalsh(effective).min() >= -1e-10


@pytest.mark.parametrize(
    "blocks",
    [
        ((0,), (0, 1)),
        ((0,),),
        ((0, 1.0), (2,)),
        ((0, 1), (2, 3)),
    ],
)
def test_collective_edge_parameter_map_rejects_non_partitions(blocks):
    from scripts.q2_relational_higher_spin_graph_qfi_gate import collective_edge_parameter_map
    with pytest.raises(ValueError, match="partition"):
        collective_edge_parameter_map(3, blocks)


def test_every_graph_reports_explicit_collective_block_pullback():
    result = run(two_j_values=(1,), trials=2)
    assert result["passed"]
    for name, row in result["graphs"]["0.5"].items():
        block = row["collective_blocking"]
        flattened = [edge for group in block["edge_partition"] for edge in group]
        assert sorted(flattened) == list(range(row["edges"]))
        assert block["effective_parameter_dimension"] == 3 * len(block["edge_partition"])
        assert block["pullback"] == "F_effective = B.T @ F_micro @ B"
        assert block["checks"]["effective_qfi_is_psd"]
        assert block["checks"]["effective_qfi_is_symmetric"]
        assert row["checks"]["collective_pullback_is_psd"]
