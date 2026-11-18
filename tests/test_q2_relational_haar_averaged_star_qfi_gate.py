import numpy as np
import pytest

from scripts.q2_relational_haar_averaged_star_qfi_gate import (
    haar_averaged_star_qfi,
    run,
    weighted_spin_ensemble,
)


@pytest.mark.parametrize("two_j", [1, 2, 3, 4])
def test_exact_haar_averaged_star_qfi(two_j):
    result = haar_averaged_star_qfi(two_j)
    assert result["passed"]
    assert result["zero_eigenvalue_count"] == 3
    assert result["positive_eigenvalue_count"] == 9
    assert result["mean_cross_link_qfi"] == pytest.approx(
        -result["expected_single_leg_qfi"] / 3.0
    )


def test_weighted_spin_ensemble_is_linear_and_psd():
    result = weighted_spin_ensemble((1, 2, 3, 4), (0.1, 0.2, 0.3, 0.4))
    assert result["passed"]
    assert result["normalized_weights"] == pytest.approx([0.1, 0.2, 0.3, 0.4])
    assert result["zero_eigenvalue_count"] == 3
    assert result["positive_eigenvalue_count"] == 9


def test_weight_normalization_is_scale_invariant():
    a = weighted_spin_ensemble((1, 2), (1.0, 3.0))
    b = weighted_spin_ensemble((1, 2), (5.0, 15.0))
    assert np.allclose(a["qfi"], b["qfi"], atol=1e-12)


@pytest.mark.parametrize("weights", [(0.0, 0.0), (-1.0, 2.0), (1.0,)])
def test_invalid_spin_weights_are_rejected(weights):
    with pytest.raises(ValueError):
        weighted_spin_ensemble((1, 2), weights)
