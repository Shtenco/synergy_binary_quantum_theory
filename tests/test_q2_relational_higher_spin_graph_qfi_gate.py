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
    assert {"path_3", "path_4", "star_4", "cycle_4"} <= set(result["graphs"]["0.5"])
    expected_boundary_legs = {
        "path_3": 8,
        "path_4": 10,
        "star_4": 10,
        "cycle_4": 8,
    }
    for spin, graph_rows in result["graphs"].items():
        for name, expected_boundary in expected_boundary_legs.items():
            row = graph_rows[name]
            assert row["boundary_half_links"] == expected_boundary
            assert row["qfi_parameter_dimension"] == 3 * row["edges"]
            assert row["qfi_rank_max_over_trials"] <= row["qfi_parameter_dimension"]
            assert row["qfi_rank_min_over_trials"] >= 0
            assert row["qfi_rank_min_over_trials"] + row["qfi_nullity_min_over_trials"] <= row["qfi_parameter_dimension"]
            assert row["checks"]["qfi_rank_nullity_partition_is_consistent"]
            assert row["qfi_spectral_maximum_over_trials"] >= row["qfi_spectral_minimum_over_trials"]

