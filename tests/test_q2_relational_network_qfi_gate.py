import numpy as np
import pytest

from scripts.q2_relational_network_qfi_gate import (
    network_relational_qfi,
    random_intertwiner,
    run,
)


def test_three_node_path_exposes_state_dependent_cross_link_qfi():
    result = run(seed=20261118, trials=24)
    assert result["passed"] is True
    assert result["gluing_probability_min"] == pytest.approx(1.0 / 16.0, abs=1e-12)
    assert result["gluing_probability_max"] == pytest.approx(1.0 / 16.0, abs=1e-12)
    assert result["diagonal_link_qfi_max_error_from_identity"] < 1e-10
    assert result["cross_block_absolute_max"] > 1e-3
    assert result["cross_block_scalar_max"] - result["cross_block_scalar_min"] > 1e-2


def test_network_qfi_is_symmetric_and_positive_semidefinite():
    rng = np.random.default_rng(1701)
    a, b, c = (random_intertwiner(rng) for _ in range(3))
    probability, qfi, residual = network_relational_qfi(a, b, c)
    assert probability == pytest.approx(1.0 / 16.0, abs=1e-12)
    assert np.allclose(qfi, qfi.T, atol=1e-11)
    assert np.linalg.eigvalsh(qfi)[0] >= -1e-10
    assert residual < 1e-10


def test_each_shared_half_link_twist_obeys_middle_node_gauss_pullthrough():
    rng = np.random.default_rng(2026)
    a, b, c = (random_intertwiner(rng) for _ in range(3))
    _, _, residual = network_relational_qfi(a, b, c)
    assert residual < 1e-10


def test_network_gate_requires_multiple_samples_for_correlation_variability():
    with pytest.raises(ValueError, match="at least two trials"):
        run(trials=1)
