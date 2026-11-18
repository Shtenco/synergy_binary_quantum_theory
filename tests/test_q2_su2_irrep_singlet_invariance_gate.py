import numpy as np
import pytest

from scripts.q2_su2_irrep_singlet_invariance_gate import (
    analyze_spin,
    invariant_singlet_tensor,
    run,
    spin_generators,
)


@pytest.mark.parametrize("two_j", [1, 2, 3, 4, 5, 6, 8])
def test_singlet_tensor_is_lie_algebra_invariant(two_j):
    tensor = invariant_singlet_tensor(two_j)
    for generator in spin_generators(two_j):
        assert np.linalg.norm(generator @ tensor + tensor @ generator.T) < 2e-12


@pytest.mark.parametrize("two_j", [1, 2, 3, 4, 5, 6, 8])
def test_singlet_tensor_has_correct_exchange_parity(two_j):
    tensor = invariant_singlet_tensor(two_j)
    assert np.allclose(tensor.T, ((-1.0) ** two_j) * tensor, atol=2e-12)


def test_finite_group_invariance_for_integer_and_half_integer_spins():
    result = run()
    assert result["passed"]
    assert result["tested_spins"] == [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]


def test_invalid_spin_label_is_rejected():
    with pytest.raises(ValueError):
        spin_generators(-1)
