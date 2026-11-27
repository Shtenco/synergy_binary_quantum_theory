import unittest
import numpy as np

from qag_postselection import (
    conditional_probe_state,
    expectation,
    source_density_from_amplitudes,
)
from qag_periodic_probe import (
    gaussian_probe_density,
    momentum_kick_unitary,
    periodic_position_momentum_operators,
    run_weak_value_kick_benchmark,
)


class QAGPeriodicProbeTests(unittest.TestCase):
    def test_momentum_operator_is_hermitian(self):
        _, momentum = periodic_position_momentum_operators(64, 20.0)
        np.testing.assert_allclose(momentum, momentum.conj().T, atol=1e-12)

    def test_kick_sign_convention(self):
        x, momentum = periodic_position_momentum_operators(128, 20.0)
        rho = gaussian_probe_density(x, sigma=1.0)
        before = expectation(momentum, rho).real
        kick = momentum_kick_unitary(x, delta_p=0.02)
        after = expectation(momentum, kick @ rho @ kick.conj().T).real
        self.assertAlmostEqual(after - before, 0.02, places=10)

    def test_weak_value_reverses_conditional_kick_in_weak_regime(self):
        result = run_weak_value_kick_benchmark(
            epsilon=0.1, delta_p_right=-1e-3
        )
        self.assertGreater(result.conditional_shift, 0.0)
        self.assertGreater(result.weak_value_prediction, 0.0)
        self.assertLess(result.relative_weak_error, 2e-4)
        self.assertAlmostEqual(
            result.success_probability,
            (np.sqrt(1.1) - np.sqrt(0.9))**2 / 4,
            places=12,
        )

    def test_decohered_source_gives_convex_branch_average(self):
        eps, delta = 0.1, -1e-3
        x, momentum = periodic_position_momentum_operators(256, 40.0)
        rho_p = gaussian_probe_density(x, sigma=1.0)
        source = np.array(
            [np.sqrt(1 + eps), np.sqrt(1 - eps)], dtype=complex
        ) / np.sqrt(2)
        rho_s = source_density_from_amplitudes(source, coherence=0.0)
        phi = np.array([1, -1], dtype=complex) / np.sqrt(2)
        ul = np.eye(len(x), dtype=complex)
        ur = momentum_kick_unitary(x, delta)
        result = conditional_probe_state(rho_s, phi, [ul, ur], rho_p)
        shift = (
            expectation(momentum, result.state)
            - expectation(momentum, rho_p)
        ).real
        self.assertLess(shift, 0.0)


if __name__ == "__main__":
    unittest.main()
