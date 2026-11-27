import unittest
import numpy as np

from qag_hamiltonian import (
    branch_unitaries_from_linear_potential,
    linear_potential_hamiltonian,
    unitary_from_diagonal_hamiltonian,
)
from qag_periodic_probe import (
    gaussian_probe_density,
    momentum_kick_unitary,
    periodic_position_momentum_operators,
)
from qag_postselection import expectation


class QAGHamiltonianTests(unittest.TestCase):
    def test_linear_potential_has_declared_force_sign(self):
        x = np.linspace(-5.0, 5.0, 32)
        h = linear_potential_hamiltonian(x, force=0.3)
        np.testing.assert_allclose(np.diag(h).real, -0.3 * x)
        np.testing.assert_allclose(h, h.conj().T, atol=1e-13)

    def test_zero_duration_is_identity_and_unitary(self):
        x = np.linspace(-5.0, 5.0, 32)
        h = linear_potential_hamiltonian(x, force=0.3)
        u = unitary_from_diagonal_hamiltonian(h, duration=0.0)
        np.testing.assert_allclose(u, np.eye(32), atol=1e-13)
        np.testing.assert_allclose(u.conj().T @ u, np.eye(32), atol=1e-13)

    def test_hamiltonian_evolution_matches_momentum_kick(self):
        n, length, force, duration = 256, 40.0, 0.02, 1.0
        x, p = periodic_position_momentum_operators(n, length)
        rho = gaussian_probe_density(x, sigma=0.8)
        u_left, u_right = branch_unitaries_from_linear_potential(
            x, force_left=0.0, force_right=force,
            duration=duration, hbar=1.0
        )
        direct = u_right @ rho @ u_right.conj().T
        kick = momentum_kick_unitary(x, force * duration)
        reference = kick @ rho @ kick.conj().T
        np.testing.assert_allclose(direct, reference, atol=1e-13)
        shift = (expectation(p, direct) - expectation(p, rho)).real
        self.assertAlmostEqual(shift, force * duration, places=10)
        np.testing.assert_allclose(u_left, np.eye(n), atol=1e-13)

    def test_rejects_non_diagonal_hamiltonian(self):
        h = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
        with self.assertRaises(ValueError):
            unitary_from_diagonal_hamiltonian(h, duration=1.0)

    def test_rejects_negative_duration(self):
        h = np.diag([0.0, 1.0]).astype(complex)
        with self.assertRaises(ValueError):
            unitary_from_diagonal_hamiltonian(h, duration=-0.1)


if __name__ == "__main__":
    unittest.main()
