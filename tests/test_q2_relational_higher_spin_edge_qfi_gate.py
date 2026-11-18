import json

import numpy as np
import pytest

from scripts.q2_relational_higher_spin_edge_qfi_gate import (
    analyze_sample,
    invariant_node_basis,
    random_invariant_node,
    run,
    spin_generators,
    _json_default,
)


@pytest.mark.parametrize("two_j", [1, 2, 3])
def test_four_valent_higher_spin_intertwiner_space_has_expected_dimension(two_j):
    basis, diagnostics = invariant_node_basis(two_j)
    assert basis.shape == ((two_j + 1) ** 4, two_j + 1)
    assert diagnostics["max_gauss_residual"] < 3e-9
    assert diagnostics["basis_orthogonality_error"] < 3e-9


@pytest.mark.parametrize("two_j", [1, 2, 3])
def test_higher_spin_two_node_qfi_matches_casimir_prediction(two_j):
    rng = np.random.default_rng(1000 + two_j)
    basis, _ = invariant_node_basis(two_j)
    node_a = random_invariant_node(rng, basis)
    node_b = random_invariant_node(rng, basis)
    result = analyze_sample(two_j, node_a, node_b, spin_generators(two_j))
    assert result["passed"]
    expected = 4 * (two_j / 2) * (two_j / 2 + 1) / 3
    assert np.asarray(result["qfi_matrix"]) == pytest.approx(expected * np.eye(3), abs=3e-9)
    assert result["gluing_probability"] == pytest.approx(1 / (two_j + 1) ** 2, abs=3e-10)


def test_higher_spin_gate_end_to_end_and_json_export():
    result = run(seed=77, trials=2)
    assert result["passed"]
    encoded = json.dumps(result, default=_json_default)
    assert '"passed": true' in encoded
