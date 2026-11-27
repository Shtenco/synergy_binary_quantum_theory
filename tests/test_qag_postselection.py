import unittest
import numpy as np

from qag_postselection import (
    complete_two_outcome_probabilities,
    conditional_probe_state,
    conditional_state_via_kraus,
    pauli_y_rotation,
    source_density_from_amplitudes,
    weak_value_right,
)


class QAGPostselectionTests(unittest.TestCase):
    def setUp(self):
        self.eps = 0.1
        self.a = np.array([
            np.sqrt(1 + self.eps), np.sqrt(1 - self.eps)
        ], dtype=complex) / np.sqrt(2)
        self.phi = np.array([1, -1], dtype=complex) / np.sqrt(2)
        self.ul = np.eye(2, dtype=complex)
        self.ur = pauli_y_rotation(0.17)
        self.rho_p = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)

    def test_weak_value_and_success_probability_reference(self):
        expected_weak = -(
            np.sqrt(1 - self.eps**2) + 1 - self.eps
        ) / (2 * self.eps)
        self.assertAlmostEqual(weak_value_right(self.a, self.phi).real,
                               expected_weak, places=12)
        # The paper's 1/eps and eps^2/4 relations are leading-order
        # asymptotics for this parameterisation, not exact identities.
        result = conditional_probe_state(
            source_density_from_amplitudes(self.a), self.phi,
            [self.ul, self.ur], self.rho_p
        )
        exact_probability = (
            np.sqrt(1 + self.eps) - np.sqrt(1 - self.eps)
        )**2 / 4
        self.assertAlmostEqual(result.probability, exact_probability, places=12)
        self.assertAlmostEqual(result.probability / (self.eps**2 / 4),
                               1.0, delta=0.06)

    def test_pure_source_kraus_matches_density_matrix_channel(self):
        rho_s = source_density_from_amplitudes(self.a, coherence=1.0)
        exact = conditional_probe_state(rho_s, self.phi,
                                        [self.ul, self.ur], self.rho_p)
        kraus = conditional_state_via_kraus(self.a, self.phi,
                                            [self.ul, self.ur], self.rho_p)
        self.assertAlmostEqual(exact.probability, kraus.probability, places=12)
        np.testing.assert_allclose(exact.state, kraus.state, atol=1e-12)

    def test_decohered_source_has_no_cross_terms(self):
        rho_s = source_density_from_amplitudes(self.a, coherence=0.0)
        exact = conditional_probe_state(rho_s, self.phi,
                                        [self.ul, self.ur], self.rho_p)
        expected = (
            abs(self.a[0])**2 * 0.5 * self.ul @ self.rho_p @ self.ul.conj().T
            + abs(self.a[1])**2 * 0.5 * self.ur @ self.rho_p @ self.ur.conj().T
        )
        np.testing.assert_allclose(exact.unnormalised, expected, atol=1e-12)

    def test_partial_coherence_interpolates_cross_terms(self):
        rho0 = source_density_from_amplitudes(self.a, coherence=0.0)
        rho1 = source_density_from_amplitudes(self.a, coherence=1.0)
        rhom = source_density_from_amplitudes(self.a, coherence=0.35)
        r0 = conditional_probe_state(rho0, self.phi, [self.ul, self.ur], self.rho_p)
        r1 = conditional_probe_state(rho1, self.phi, [self.ul, self.ur], self.rho_p)
        rm = conditional_probe_state(rhom, self.phi, [self.ul, self.ur], self.rho_p)
        # Unnormalised channel is affine in the source coherence.
        np.testing.assert_allclose(
            rm.unnormalised,
            r0.unnormalised + 0.35 * (r1.unnormalised-r0.unnormalised),
            atol=1e-12
        )

    def test_two_outcome_probabilities_sum_to_one(self):
        rho_s = source_density_from_amplitudes(self.a, coherence=0.42)
        p0, p1 = complete_two_outcome_probabilities(
            rho_s, self.rho_p, [self.ul, self.ur]
        )
        self.assertAlmostEqual(p0 + p1, 1.0, places=12)

    def test_normalisation_and_positivity(self):
        rho_s = source_density_from_amplitudes(self.a, coherence=0.7)
        result = conditional_probe_state(rho_s, self.phi,
                                         [self.ul, self.ur], self.rho_p)
        self.assertAlmostEqual(np.trace(result.state).real, 1.0, places=12)
        self.assertGreaterEqual(np.linalg.eigvalsh(result.state).min(), -1e-12)


if __name__ == "__main__":
    unittest.main()
