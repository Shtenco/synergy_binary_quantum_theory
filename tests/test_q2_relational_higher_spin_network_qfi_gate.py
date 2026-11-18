import numpy as np
import pytest

from scripts.q2_relational_higher_spin_network_qfi_gate import (
    analyze_spin,
    four_valent_intertwiner_basis,
    invariant_singlet_tensor,
    run,
    spin_generators,
)


@pytest.mark.parametrize("two_j", [1, 2, 3, 4])
def test_four_valent_intertwiner_dimension_and_gauss_law(two_j):
    basis, residual = four_valent_intertwiner_basis(two_j)
    assert basis.shape == ((two_j + 1) ** 4, two_j + 1)
    assert residual < 2e-10


@pytest.mark.parametrize("two_j", [1, 2, 3, 4])
def test_spin_j_singlet_invariance(two_j):
    singlet = invariant_singlet_tensor(two_j)
    for generator in spin_generators(two_j):
        assert np.linalg.norm(generator @ singlet + singlet @ generator.T) < 2e-10


@pytest.mark.parametrize("two_j", [1, 2, 3])
def test_higher_spin_path_qfi_identities(two_j):
    result = analyze_spin(two_j, seed=20261118, trials=3)
    assert all(result["checks"].values())


def test_higher_spin_gate_runs_half_integer_and_integer_spins():
    result = run(two_j_values=(1, 2, 3), trials=2)
    assert result["passed"]
    assert result["tested_spins"] == [0.5, 1.0, 1.5]


def test_invalid_spin_label_is_rejected():
    with pytest.raises(ValueError):
        spin_generators(0)
