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
