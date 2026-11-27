import unittest
import numpy as np

from qag_relational_apparatus import (
    apply_relational_interaction,
    cyclic_translation,
    gaussian_position_state,
    relative_interaction_unitary,
    translation_commutator_norm,
)


class QAGRelationalApparatusTests(unittest.TestCase):
    def test_relative_unitary_is_unitary(self):
        unitary = relative_interaction_unitary(12, coupling=0.7, duration=1.3)
        np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(144), atol=1e-12)

    def test_simultaneous_translation_commutes_exactly(self):
        unitary = relative_interaction_unitary(10, coupling=1.2, duration=0.8, harmonic=2)
        self.assertLess(translation_commutator_norm(unitary, 10), 1e-12)

    def test_translation_operator_is_unitary(self):
        shift = cyclic_translation(9, steps=3)
        np.testing.assert_allclose(shift.conj().T @ shift, np.eye(9), atol=1e-12)

    def test_interaction_preserves_norm_and_translation_expectation(self):
        source = gaussian_position_state(16, center=3.0, width=1.2)
        apparatus = gaussian_position_state(16, center=11.0, width=1.5)
        result = apply_relational_interaction(source, apparatus, 0.9, 1.1)
        self.assertAlmostEqual(result.norm, 1.0, places=12)
        self.assertAlmostEqual(result.translation_expectation_initial.real,
                               result.translation_expectation_final.real, places=12)
        self.assertAlmostEqual(result.translation_expectation_initial.imag,
                               result.translation_expectation_final.imag, places=12)
        self.assertLess(result.commutator_norm, 1e-12)

    def test_nonzero_relational_coupling_generates_entanglement(self):
        source = gaussian_position_state(16, center=3.0, width=2.0)
        apparatus = gaussian_position_state(16, center=11.0, width=2.0)
        uncoupled = apply_relational_interaction(source, apparatus, 0.0, 1.0)
        coupled = apply_relational_interaction(source, apparatus, 1.7, 1.0)
        self.assertLess(uncoupled.entanglement_entropy_bits, 1e-10)
        self.assertGreater(coupled.entanglement_entropy_bits, 1e-5)
        self.assertAlmostEqual(np.trace(coupled.source_density).real, 1.0, places=12)
        self.assertAlmostEqual(np.trace(coupled.apparatus_density).real, 1.0, places=12)

    def test_joint_readout_probabilities_are_normalized(self):
        source = gaussian_position_state(12, center=2.0, width=1.0)
        apparatus = gaussian_position_state(12, center=8.0, width=1.0)
        result = apply_relational_interaction(source, apparatus, 0.8, 0.9)
        self.assertAlmostEqual(float(result.joint_position_probability.sum()), 1.0, places=12)
        self.assertTrue(np.all(result.joint_position_probability >= 0.0))

    def test_rejects_invalid_parameters(self):
        with self.assertRaises(ValueError):
            gaussian_position_state(3, 0.0, 1.0)
        with self.assertRaises(ValueError):
            relative_interaction_unitary(8, 1.0, 1.0, hbar=0.0)
        with self.assertRaises(ValueError):
            relative_interaction_unitary(8, 1.0, 1.0, harmonic=0)


if __name__ == "__main__":
    unittest.main()
