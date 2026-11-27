import unittest
import numpy as np
from qag_measurement_dilation import (
    apply_source_measurement_dilation,
    measurement_outcome_probabilities,
    source_measurement_dilation_unitary,
)

class QAGMeasurementDilationTests(unittest.TestCase):
    def test_measurement_unitary_is_exactly_unitary(self):
        theta = 0.73
        basis = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]], dtype=complex)
        unitary = source_measurement_dilation_unitary(basis)
        np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(4), atol=1e-12)

    def test_meter_probabilities_match_born_rule(self):
        basis = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
        source = np.array([np.sqrt(0.8), np.sqrt(0.2)], dtype=complex)
        dilated = apply_source_measurement_dilation(source, basis)
        actual = measurement_outcome_probabilities(dilated)
        expected = tuple(abs(np.vdot(basis[:, k], source)) ** 2 for k in range(2))
        np.testing.assert_allclose(actual, expected, atol=1e-12)
        self.assertAlmostEqual(sum(actual), 1.0, places=12)

    def test_entangled_environment_is_retained_until_meter_readout(self):
        basis = np.eye(2, dtype=complex)
        state = np.array([[1, 0], [0, 1]], dtype=complex) / np.sqrt(2)
        dilated = apply_source_measurement_dilation(state, basis)
        self.assertEqual(dilated.shape, (2, 2, 2))
        p0, p1 = measurement_outcome_probabilities(dilated)
        self.assertAlmostEqual(p0, 0.5, places=12)
        self.assertAlmostEqual(p1, 0.5, places=12)
        np.testing.assert_allclose(dilated[:, :, 0], np.diag([1, 0]) / np.sqrt(2))
        np.testing.assert_allclose(dilated[:, :, 1], np.diag([0, 1]) / np.sqrt(2))

    def test_meter_register_is_not_silently_discarded(self):
        source = np.array([1, 1], dtype=complex) / np.sqrt(2)
        dilated = apply_source_measurement_dilation(source, np.eye(2))
        self.assertAlmostEqual(np.linalg.norm(dilated), 1.0, places=12)
        self.assertAlmostEqual(np.vdot(dilated[..., 0], dilated[..., 1]), 0.0, places=12)

    def test_rejects_invalid_basis_and_state(self):
        with self.assertRaises(ValueError): source_measurement_dilation_unitary(np.ones((2, 2)))
        with self.assertRaises(ValueError): apply_source_measurement_dilation(np.ones(3), np.eye(2))
        with self.assertRaises(ValueError): apply_source_measurement_dilation(np.array([1, np.inf]), np.eye(2))

if __name__ == "__main__":
    unittest.main()