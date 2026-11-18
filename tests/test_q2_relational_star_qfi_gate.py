import numpy as np
import pytest

from scripts.q2_relational_star_qfi_gate import pair_casimir_matrix, random_intertwiner, run


def test_four_valent_star_qfi_has_three_gauss_gauge_null_directions():
    result = run(seed=20261118, trials=24)
    assert result["passed"] is True
    assert result["local_parameter_count"] == 12
    assert result["max_common_twist_null_residual"] < 1e-11
    assert result["max_full_qfi_gauge_null_residual"] < 1e-10
    assert result["minimum_gauge_fixed_relative_metric_eigenvalue"] > 1e-10


def test_pair_metric_is_symmetric_psd_and_has_unit_diagonal():
    rng = np.random.default_rng(991)
    for _ in range(10):
        state = random_intertwiner(rng)
        metric = pair_casimir_matrix(state)
        assert np.allclose(metric, metric.T, atol=1e-11)
        assert np.allclose(np.diag(metric), np.ones(4), atol=1e-11)
        assert np.linalg.eigvalsh(metric)[0] >= -1e-10
        assert np.allclose(metric @ np.ones(4), np.zeros(4), atol=1e-11)


def test_fixing_one_link_removes_the_common_twist_null_mode():
    rng = np.random.default_rng(1204)
    state = random_intertwiner(rng)
    metric = pair_casimir_matrix(state)
    relative = metric[:3, :3]
    assert np.linalg.eigvalsh(relative)[0] > 1e-10


def test_star_gate_requires_multiple_samples():
    with pytest.raises(ValueError, match="at least two trials"):
        run(trials=1)
