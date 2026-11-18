import numpy as np
import pytest

from scripts.q2_relational_edge_qfi_gate import (
    random_intertwiner,
    relational_boundary_state,
    relational_qfi_matrix,
    run,
)


def test_two_block_singlet_gluing_has_unit_relational_qfi():
    result = run(seed=20261118, trials=12)
    assert result["passed"] is True
    assert result["gluing_probability_min"] == pytest.approx(0.25, abs=1e-12)
    assert result["gluing_probability_max"] == pytest.approx(0.25, abs=1e-12)
    assert np.allclose(result["relational_qfi_mean"], np.eye(3), atol=1e-12)
    assert np.allclose(result["normalized_qfi_metric_F_over_4"],
                       0.25 * np.eye(3), atol=1e-12)


def test_relational_qfi_is_independent_of_the_two_node_intertwiners():
    rng = np.random.default_rng(71)
    for _ in range(8):
        a = random_intertwiner(rng)
        b = random_intertwiner(rng)
        probability, qfi, pullthrough_residual = relational_qfi_matrix(a, b)
        assert probability == pytest.approx(0.25, abs=1e-12)
        assert np.allclose(qfi, np.eye(3), atol=1e-11)
        assert pullthrough_residual < 1e-11


def test_relative_twist_is_not_the_global_node_gauss_generator():
    # Each input block is a four-valent singlet, but the inserted generator
    # acts on one shared half-link only (not on the sum over all four legs).
    result = run(seed=19, trials=2)
    assert result["max_node_gauss_residual"] < 1e-12
    assert result["max_half_link_to_boundary_generator_residual"] < 1e-11
    assert result["relational_qfi_max_frobenius_error_from_identity"] < 1e-11

def test_relational_density_operator_is_normalized_and_pure():
    rng = np.random.default_rng(90210)
    a = random_intertwiner(rng)
    b = random_intertwiner(rng)
    state, probability = relational_boundary_state(a, b)
    rho_e = np.outer(state, state.conj())
    assert probability == pytest.approx(0.25, abs=1e-12)
    assert np.trace(rho_e) == pytest.approx(1.0 + 0.0j, abs=1e-12)
    assert np.trace(rho_e @ rho_e) == pytest.approx(1.0 + 0.0j, abs=1e-12)
